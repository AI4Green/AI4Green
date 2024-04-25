""" For tabbular data export formats - CSV or SURF"""
import io
import json
import math
from typing import Dict, List, Literal, Optional

import azure.core.exceptions
import pandas as pd
from flask import abort
from sources import models, services


class CsvExport:
    """CSV export intended to have all the data to recreate an experiment in AI4Green"""

    def __init__(self, data_export_request: models.DataExportRequest):
        self.data_export_request = data_export_request

    def make_row(self, reaction):
        # first get the data included in the surf export
        data = SurfExport(self.data_export_request).make_row(reaction).to_dict()
        # now get the data not included within the surf export
        data["time_of_creation"] = reaction.time_of_creation
        data["time_of_update"] = reaction.time_of_update
        data["complete"] = reaction.complete
        # get the sustainability data
        metadata_dict = services.data_export.metadata.ReactionMetaData(
            reaction
        ).make_metadata_dict()
        # combine these dicts without duplicating entries.
        print(metadata_dict)


class SurfExport:
    amount_factor_dict = {"μmol": 0.000001, "mmol": 0.001, "mol": 1}

    """
    Format described in more detail here https://github.com/alexarnimueller/surf.
    Functionality to convert SURF to open reaction database format can also be found at the SURF github.
    """

    def __init__(self, data_export_request: models.DataExportRequest):
        self.data_export_request = data_export_request

    def make_row(self, reaction: models.Reaction):
        summary_data = json.loads(reaction.summary_table_data)
        reaction_data = json.loads(reaction.reaction_table_data)
        # we need to get all the necessary data and add it to a dict
        data = {}
        data["rxn_id"] = reaction.reaction_id
        data["source_id"] = "AI4Green.app"
        data["source_type"] = "ELN"
        data["rxn_type"] = None
        data["temperature_deg_c"] = self._read_value(
            summary_data, "reaction_temperature"
        )
        data["time_h"] = None
        data["atmosphere"] = None
        data["stirring_shaking"] = None
        data["scale_mol"] = self._get_mol_scale(reaction_data)
        data["concentration_mol_l"] = self._get_concentration(
            reaction_data.get("solvent_concentrations")
        )
        self._add_reagent_and_catalyst_data(data, reaction, reaction_data)
        self._add_solvent_data(data, reaction, reaction_data)
        self._add_starting_material_data(data, reaction, reaction_data)
        self._add_product_data(data, reaction, reaction_data, summary_data)
        creator = reaction.creator_person.user
        data[
            "provenance"
        ] = f"{creator.fullname} {creator.email} {reaction.workbook.WorkGroup.name}"
        data["procedure"] = reaction.description
        df = pd.DataFrame(data, index=[0])
        return df

    def _add_product_data(
        self,
        data: Dict,
        reaction: models.Reaction,
        reaction_data: Dict,
        summary_data: Dict,
    ):
        """
        Updates the data export dict with product data
        Args:
            data - the data export dictionary which is updated in the function
            reaction - the reaction we are getting data for
            reaction_data - the reaction table data for the reaction we are looking at.
        """
        if services.data_export.utils.check_compounds_present(reaction, "products"):
            for idx, product_smiles in enumerate(reaction.products, 1):
                data[f"product_{idx}_cas"] = services.all_compounds.cas_from_smiles(
                    product_smiles
                )
                data[f"product_{idx}_ms"] = None
                data[f"product_{idx}_nmr"] = None
                data[f"product_{idx}_smiles"] = product_smiles
                # only main product yield is recorded
                if str(idx) == reaction_data["main_product"]:
                    percent_yield = self._get_percent_yield(reaction_data, summary_data)
                else:
                    percent_yield = None
                data[f"product_{idx}_yield"] = percent_yield
                data[f"product_{idx}_yieldtype"] = None

    @staticmethod
    def _add_starting_material_data(
        data: Dict, reaction: models.Reaction, reaction_data: Dict
    ):
        """
        Updates the data export dict with reactant/starting material data
        Args:
            data - the data export dictionary which is updated in the function
            reaction - the reaction we are getting data for
            reaction_data - the reaction table data for the reaction we are looking at.
        """
        if services.data_export.utils.check_compounds_present(reaction, "reactants"):
            for idx, reactant_smiles in enumerate(reaction.reactants, 1):
                data[f"startingmat_{idx}_eq"] = reaction_data["reactant_equivalents"][
                    idx - 1
                ]
                data[f"startingmat_{idx}_cas"] = services.all_compounds.cas_from_smiles(
                    reactant_smiles
                )
                data[f"startingmat_{idx}_smiles"] = reactant_smiles

    @staticmethod
    def _add_solvent_data(data: Dict, reaction: models.Reaction, reaction_data: Dict):
        """
        Updates the data export dict with solvent data
        Args:
            data - the data export dictionary which is updated in the function
            reaction - the reaction we are getting data for
            reaction_data - the reaction table data for the reaction we are looking at.
        """
        if services.data_export.utils.check_compounds_present(reaction, "solvent"):
            for idx, solvent_pk in enumerate(reaction.solvent, 1):
                data[f"solvent_{idx}_cas"] = services.solvent.cas_from_primary_key(
                    solvent_pk
                )
                # get volumes of all solvents in this reaction
                total_solvent_volume = sum(
                    [float(x) for x in reaction_data["solvent_volumes"]]
                )
                data[f"solvent_{idx}_fraction"] = (
                    float(reaction_data["solvent_volumes"][idx - 1])
                    / total_solvent_volume
                )
                data[
                    f"solvent_{idx}_smiles"
                ] = services.solvent.smiles_from_primary_key(solvent_pk)

    @staticmethod
    def _get_percent_yield(reaction_data: Dict, summary_data: Dict) -> Optional[float]:
        actual_mass_yield = summary_data.get("real_product_mass")
        # a yield of zero is a valid yield
        if actual_mass_yield:
            percent_yield = services.sustainability.percent_yield_from_rxn_data(
                reaction_data, float(actual_mass_yield)
            )
        else:
            percent_yield = None
        return percent_yield

    @staticmethod
    def _reagent_or_catalyst(equivalent: float) -> Literal["reagent", "catalyst"]:
        """
        Uses the equivalent of a reagent to determine to label it a reagent or a catalyst

        """
        if equivalent > 0.25:
            return "reagent"
        return "catalyst"

    @staticmethod
    def _get_concentration(solvent_concentrations: List) -> Optional[float]:
        """Gets the concentration the reaction was performed at using volume of solvent and limited reactants."""
        if not solvent_concentrations:
            return
        solvent_concentrations = [float(x) for x in solvent_concentrations]
        return math.prod(solvent_concentrations)

    def _get_mol_scale(self, reaction_data: Dict) -> Optional[float]:
        """
        To get the mol scale, look up the limiting reactant and get this amount and its units and convert to mol
        Args:
            reaction_data - the reaction table data for the reaction we are looking at.
        Returns:
            The mol scale of the reaction expressed as a float
        """
        limiting_reactant = reaction_data.get("limiting_reactant_table_number")
        if limiting_reactant and reaction_data.get("reactant_amounts_raw"):
            limiting_reactant_idx = int(limiting_reactant) - 1
            limiting_reactant_amount = round(
                float(reaction_data.get("reactant_amounts_raw")[limiting_reactant_idx]),
                2,
            )
            # convert to mol scale
            reactant_amount_unit = reaction_data["amount_units"]
            amount_factor = self.amount_factor_dict[reactant_amount_unit]
            return limiting_reactant_amount * amount_factor
        return None

    @staticmethod
    def _read_value(summary_data, key) -> Optional[str]:
        """Returns the reaction temperature the user entered into the summary table"""
        value = summary_data.get(key)
        return value if value != "" else None

    def _add_reagent_and_catalyst_data(
        self, data: Dict, reaction: models.Reaction, reaction_data: Dict
    ):
        """
        Updates the data export dict with catalyst and reagent data
        Args:
            data - the data export dictionary which is updated in the function
            reaction - the reaction we are getting data for
            reaction_data - the reaction table data for the reaction we are looking at.
        """
        if services.data_export.utils.check_compounds_present(reaction, "reagents"):
            reagent_idx_list = []
            catalyst_idx_list = []
            for idx, reagent_equivalent in enumerate(
                reaction_data.get("reagent_equivalents")
            ):
                role = self._reagent_or_catalyst(float(reagent_equivalent))
                if role == "reagent":
                    reagent_idx_list.append(idx)
                if role == "catalyst":
                    catalyst_idx_list.append(idx)
            self._update_reagent_or_catalyst_dict(
                data, reaction, reaction_data, "reagent", reagent_idx_list
            )
            self._update_reagent_or_catalyst_dict(
                data, reaction, reaction_data, "catalyst", catalyst_idx_list
            )

    @staticmethod
    def _update_reagent_or_catalyst_dict(
        data: Dict,
        reaction: models.Reaction,
        reaction_data: Dict,
        role: str,
        compound_idx_list: List[int],
    ):
        if compound_idx_list:
            for label_idx, idx in enumerate(compound_idx_list, 1):
                compound_equivalent = reaction_data["reagent_equivalents"][idx]
                compound_smiles = reaction.reagents[idx]
                data[f"{role}_{label_idx}_eq"] = compound_equivalent
                data[
                    f"{role}_{label_idx}_cas"
                ] = services.all_compounds.cas_from_smiles(compound_smiles)
                data[f"{role}_{label_idx}_smiles"] = compound_smiles

    def save(self, df: pd.DataFrame, filename: str):
        """Calls the appropriate method to save the data. Can't use .RDF without a reaction_container object"""
        """Saves the blob to Azure blob service"""

        csv_buffer = io.StringIO()
        df.to_csv(csv_buffer, index=False, sep="\t")
        # Get CSV file contents from csv buffer and save as a csv in the export container in Azure
        file_contents = bytearray(csv_buffer.getvalue(), "utf-8")

        # Upload the blob data
        upload = io.BytesIO(file_contents)
        blob_client = services.file_attachments.get_blob_client(
            "export-outputs", filename
        )

        # todo remove before deployment because this should only happen when debugging
        blob_client.upload_blob(upload, blob_type="BlockBlob", overwrite=True)

        # confirm upload
        if not blob_client.exists():
            print(f"blob {filename} upload failed")
            abort(401)
