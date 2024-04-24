""" For tabbular data export formats - CSV or SURF"""
import json
import math
from typing import Dict, List, Literal, Optional

from sources import models, services


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
        data["temperature_deg_c"] = self._read_value(summary_data, "temperature")
        data["time_h"] = None
        data["atmosphere"] = None
        data["stirring_shaking"] = None
        data["scale_mol"] = self._get_mol_scale(reaction_data)
        data["concentration_mol_l"] = self._get_concentration(
            reaction_data.get("solvent_concentrations")
        )
        self._add_reagent_and_catalyst_data(data, reaction, reaction_data)
        if services.data_export.utils.check_compounds_present(reaction, "solvent"):
            for idx, solvent in enumerate(reaction.solvent, 1):
                solvent_db = services.solvent.from_primary_key(solvent)
                print(solvent_db)
                # cas = solvent_db.cas if solvent_db else None
                # data[f"solvent_{idx}_cas"] = solvent_db.cas

                # solvent_volume = reaction_data

                # get data from primary key

        # reagent and catalyst are same
        # data['catalyst_1_cas']
        # data['catalyst_1_eq']
        # data['catalyst_1_smiles']
        #
        # # solvent
        # ['cas', 'fraction', 'SMILES']
        #
        #
        # # startingmat
        # ['cas', 'eq', 'smiles']

    def _reagent_or_catalyst(self, equivalent: float) -> Literal["reagent", "catalyst"]:
        """
        Uses the equivalent of a reagent to determine to label it a reagent or a catalyst

        """
        if equivalent > 0.25:
            return "reagent"
        return "catalyst"

    def _get_concentration(self, solvent_concentrations: List) -> Optional[float]:
        """Gets the concentration the reaction was performed at using volume of solvent and limited reactants."""
        if (
            not isinstance(solvent_concentrations, list)
            or not solvent_concentrations[0].isdigit()
        ):
            return
        solvent_concentrations = [float(x) for x in solvent_concentrations]
        return math.prod(solvent_concentrations)

    def _get_mol_scale(self, reaction_data: Dict) -> Optional[float]:
        limiting_reactant = reaction_data.get("limiting_reactant_table_number")
        if limiting_reactant:
            limiting_reactant_idx = int(limiting_reactant) - 1
            if reaction_data.get("reactant_amounts_raw"):
                limiting_reactant_amount = round(
                    float(
                        reaction_data.get("reactant_amounts_raw")[limiting_reactant_idx]
                    ),
                    2,
                )
                # convert to mol scale
                reactant_amount_unit = reaction_data["amount_units"]
                amount_factor = self.amount_factor_dict[reactant_amount_unit]
                return limiting_reactant_amount * amount_factor
        return None

    def _read_value(self, summary_data, key) -> Optional[str]:
        """Returns the reaction temperature the user entered into the summary table"""
        value = summary_data.get("reaction_temperature")
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

            if catalyst_idx_list:
                for catalyst_idx, idx in enumerate(catalyst_idx_list, 1):
                    catalyst_equivalent = reaction.reagent_equivalents[idx]
                    catalyst_smiles = reaction.reagents[idx]
                    data[f"catalyst_{catalyst_idx}_eq"] = catalyst_equivalent
                    data[
                        f"catalyst_{catalyst_idx}_cas"
                    ] = services.all_compounds.cas_from_smiles(catalyst_smiles)
                    data[f"catalyst_{catalyst_idx}_smiles"] = catalyst_smiles

            if reagent_idx_list:
                for reagent_idx, idx in enumerate(reagent_idx_list, 1):
                    reagent_equivalent = reaction.reagent_equivalents[idx]
                    reagent_smiles = reaction.reagents[idx]
                    data[f"reagent_{reagent_idx}_eq"] = reagent_equivalent
                    data[
                        f"reagent_{reagent_idx}_cas"
                    ] = services.all_compounds.cas_from_smiles(reagent_smiles)
                    data[f"reagent_{reagent_idx}_smiles"] = reagent_smiles


# def _save_blob(self):
#     """Saves the blob to Azure blob service"""
#     blob_client = services.file_attachments.get_blob_client(
#         self.container_name, self.filename
#     )
#     # Upload the blob data
#     upload = io.BytesIO(self.file_contents)
#
#     try:  # todo remove before deployment because this should only happen when debugging
#         blob_client.upload_blob(upload, blob_type="BlockBlob")
#     except azure.core.exceptions.ResourceExistsError:
#         print("duplicate resource name error. skipping second one")
#         pass
#
#     # confirm upload
#     if not blob_client.exists():
#         print(f"blob {self.filename} upload failed")
#         abort(401)
