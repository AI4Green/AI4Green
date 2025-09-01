import re
from typing import List, Tuple

from flask import json, jsonify, render_template
from sources import services
from sources.extensions import db


def get_reactants_and_products_list(
    reaction_smiles: str,
) -> Tuple[List[str], List[str]]:
    """
    Process reactants and products strings to obtain lists of reactant and product SMILES.
    Converts words into SMILES symbols and identifies the format as either CXSMILES or SMILES.
    If ions are present, it processes the ionic SMILES accordingly.

    Args:
        reaction_smiles (str), the smiles of the reaction from the front end

    Returns:
        Tuple[List[str], List[str]]: A tuple containing lists of reactant, reagent and product SMILES.

    """
    # reactants_smiles = smiles_symbols(reactants)
    # products_smiles = smiles_symbols(products)

    # # form the reaction_smiles. Just from reactands and products to exclude reagents/other data
    # reaction_smiles = (reactants_smiles + ">>" + products_smiles).replace(",", ".")

    # we get the smiles straight from the sketcher to see if we have CXSmiles
    # sketcher_smiles = smiles_symbols(request.json.get("reaction_smiles"))

    # [OH-].[Na+]>>[Cl-].[Cl-].[Zn++] |f:0.1,2.3.4|
    if "|f:" in reaction_smiles:
        reaction_smiles += re.search(r" \|[^\|]*\|$", reaction_smiles).group()
        (
            reactants_smiles_list,
            products_smiles_list,
        ) = services.ions.reactants_and_products_from_ionic_cx_smiles(reaction_smiles)
    elif re.findall(r"(?<!\{)\+", reaction_smiles) or re.findall(
        r"(?<!\{)-", reaction_smiles
    ):
        # find "+" or "-" in reaction_smiles but not {-} or {+n} for polymers
        # CHANGE THIS FUNCTION TO HANDLE REAGENTS
        (
            reactants_smiles_list,
            products_smiles_list,
        ) = services.ions.reactants_and_products_from_ionic_smiles(reaction_smiles)
    # reactions with no ions - make rxn object directly from string
    else:
        # Split and pad to ensure two parts
        # change to 3 for reagent support and change split type >> to >
        parts = reaction_smiles.split(">>") + [""] * (
            2 - len(reaction_smiles.split(">>"))
        )

        # Process each part (add reagent_smiles_list_here)
        reactants_smiles_list, products_smiles_list = [
            x.split(".") if x else [] for x in parts
        ]

    return reactants_smiles_list, products_smiles_list


class ReactionTable:
    """
    Class that stores all the associated data for a reaction table
    Made to update reaction table row by row instead of erasing all data
    not sure if this will work, but yolo lmao
    """

    def __init__(self, reaction, workgroup, workbook, demo, tutorial):
        # set up default values
        self.reaction = reaction
        self.reaction_table_data = {}
        self.units = self.default_units()
        self.workgroup = workgroup
        self.workbook = workbook
        self.demo = demo
        self.tutorial = tutorial

        # to do
        # fix polymer bugs
        # fix input novel compound

        # reaction table species
        self.reactants = []
        self.products = []
        self.reagents = []
        self.solvents = []

        # interactive elements
        self.solvent_dropdown = self.get_solvent_dropdown()
        self.reaction_class = None

        # novel_compound responses are stored here
        self.novel_compounds = []

        # params/workups/purifications/responses
        self.reaction_params = {}

        # load reaction_table from provided reaction_table data
        self.load()

    def load(self):
        self.reaction_table_data = json.loads(self.reaction.reaction_table_data)
        compounds = services.compound.SketcherCompound.from_reaction_table_dict(
            self.reaction_table_data, self.workbook
        )
        self.reactants = compounds["reactant"]
        self.products = compounds["product"]
        self.reagents = compounds["reagent"]
        self.solvents = compounds["solvent"]

        # TODO: we should think about adding new relational db tables for workups/purifications/responses
        # TODO: not sure if this would be worth it, but would add more structure

        self.load_units()
        self._load_reaction_params()
        self._load_responses()

    def check_errors(self):
        for compound_group in [self.reactants, self.products]:
            errors = self._check_compound_errors(compound_group)
            if errors:
                return errors
        return False

    @staticmethod
    def _check_compound_errors(compound_type):
        for compound in compound_type:
            if compound.novel_compound_table:
                return compound.novel_compound_table
            elif compound.errors:
                return compound.errors[0]
        return False

    def load_units(self):
        # response units not included here, as these are inputted by users
        unit_keys = [
            "limiting_reactant_table_number",
            "main_product",
            "mass_units",
            "polymerisation_type",
            "amount_units",
            "volume_units",
            "solvent_volume_units",
            "product_mass_units",
            "product_amount_units",
            "reaction_temperature_units",
            "reaction_time_units",
        ]

        self.units.update(
            {
                key: self.reaction_table_data[key]
                for key in unit_keys
                if key in self.reaction_table_data
            }
        )

    def update(self, reaction_smiles, polymer_indices):
        (
            reactants_smiles_list,
            product_smiles_list,
        ) = services.reaction_table.get_reactants_and_products_list(reaction_smiles)

        # only classify reaction if not polymer
        if not polymer_indices:
            self.reaction_class = services.reaction_classification.classify_reaction(
                reactants_smiles_list, product_smiles_list
            )

        updated_reactants = self.convert_smiles_to_compound(
            reactants_smiles_list, "Reactant", polymer_indices, 0
        )

        updated_products = self.convert_smiles_to_compound(
            product_smiles_list, "Product", polymer_indices, len(updated_reactants)
        )

        self.update_reactants_and_products(updated_reactants, updated_products)
        self.update_reaction_table_data()
        self.check_errors()

    def update_reaction_table_data(self):
        """
        updates reaction_table_data dict with info from reactionTable class
        """
        self.update_dict_reaction_component(self.reactants)
        self.update_dict_reaction_component(self.products)
        self.update_dict_reaction_component(self.reagents)
        self.update_dict_reaction_component(self.solvents)

    def save(self):
        """
        Updates the reaction_table_data dict and saves to the reaction.
        used only for reactwise imports, as reaction autosave manually saves reaction_table_dict from user input
        """
        self.reaction.update(**{"reaction_table_data": self.reaction_table_data})
        db.session.commit()

    def update_dict_reaction_component(self, reaction_component):
        """
        Updates self.reaction_table_data with names from provided reaction component.
        """
        # define which fields you want to collect
        fields = {
            "names": [],
            "smiles": [],
            "molecular_weights": [],
            "hazards": [],
            "density": [],
            "equivalents": [],
        }

        prefix = ""

        # collect values into lists
        for c in reaction_component:
            # figure out component prefix (e.g. "reactants", "products", ...)
            prefix = c.reaction_component.lower()
            # append smiles to novel compounds if novel compound
            if c.is_novel_compound:
                self.novel_compounds.append(c.smiles)
            for key, acc in fields.items():
                # should maintain smiles even if novel_compound
                if key == "smiles":
                    acc.append(c.smiles)
                else:
                    acc.extend(c.compound_data.get(key, [""]))

        # update reaction_table_data in one go
        for field, values in fields.items():
            if values:
                self.reaction_table_data[f"{prefix}_{field}"] = values

    def update_reactants_and_products(self, updated_reactants, updated_products):
        new_reactants = self.identify_compounds(updated_reactants, self.reactants)
        new_products = self.identify_compounds(updated_products, self.products)

        self.reactants = new_reactants
        self.products = new_products

    def identify_compounds(self, updated_compounds, original_compounds):
        added_compounds = self.find_added_compounds(
            updated_compounds, original_compounds
        )
        preserved_compounds = self.find_preserved_compounds(
            updated_compounds, original_compounds
        )
        return preserved_compounds + added_compounds

    def convert_smiles_to_compound(
        self, smiles_list, reaction_component, polymer_indices, number_of_reactants
    ):
        base_idx = 1 + number_of_reactants

        return [
            services.compound.SketcherCompound(
                smiles=smiles,
                idx=base_idx + i,
                polymer_indices=polymer_indices,
                workbook=self.workbook,
                demo=self.demo,
                reaction_component=reaction_component,
                reaction_component_idx=i + 1,
            )
            for i, smiles in enumerate(smiles_list)
        ]

    def _load_reaction_params(self):
        param_keys = [
            "reaction_time",
            "reaction_atmosphere",
            "reaction_temperature",
        ]
        for param in param_keys:
            # important to add try loop for older reactions
            try:
                self.reaction_params[param] = self.reaction_table_data[param]
            except KeyError:
                self.reaction_params[param] = ""

    def _load_responses(self):
        # TODO: If responses are added to the db they should be loaded in this fn
        response_names = self.reaction_table_data.get("response_names", [])
        response_units = self.reaction_table_data.get("response_units", [])
        response_values = self.reaction_table_data.get("response_values", [])
        response_notes = self.reaction_table_data.get("response_notes", [])

        self.responses = [
            {"name": n, "unit": u, "value": v, "note": nt}
            for n, u, v, nt in zip(
                response_names, response_units, response_values, response_notes
            )
        ]

    def render(self):
        reaction_table = render_template(
            "reactions/_reaction_table.html",
            reactants=self.reactants,
            reagents=self.reagents,
            solvents=self.solvents,
            products=self.products,
            reaction_params=self.reaction_params,
            responses=self.responses,
            number_of_reactants=len(self.reactants),
            number_of_reagents=len(self.reagents),
            number_of_solvents=len(self.solvents),
            number_of_products=len(self.products),
            number_of_responses=len(self.responses),
            units=self.units,
            # do we still need these?
            identifiers=[],
            reactant_table_numbers=[],
            # product_intended_dps=product_data["intended_dps"],
            reagent_table_numbers=[],
            reaction_table_data="",
            summary_table_data="",
            sol_rows=self.solvent_dropdown,
            reaction=self.reaction,
            # need to update reaction class as class attr?
            reaction_class=self.reaction.reaction_class,
            # polymer indices?
            polymer_indices={},
        )
        return jsonify({"reactionTable": reaction_table})

    @staticmethod
    def default_units():
        return {
            "limiting_reactant_table_number": 1,
            "main_product": -1,
            "mass_units": "mg",
            "polymerisation_type": "NAN",
            "amount_units": "mmol",
            "volume_units": "mL",
            "solvent_volume_units": "mL",
            "product_mass_units": "mg",
            "product_amount_units": "mmol",
            "reaction_temperature_units": "°C",
            "reaction_time_units": "h",
        }

    @staticmethod
    def find_added_compounds(updated_compounds, original_compounds):
        """
        returns unique compounds from updated compound list by comparing to current reaction table data
        """
        original_inchis = {x.inchi for x in original_compounds}
        return [x for x in updated_compounds if x.inchi not in original_inchis]

    @staticmethod
    def find_preserved_compounds(updated_compounds, original_compounds):
        updated_inchis = {x.inchi for x in updated_compounds}
        return [x for x in original_compounds if x.inchi in updated_inchis]

    def get_solvent_dropdown(self):
        if self.demo == "demo":
            sol_rows = services.solvent.get_default_list()
        else:
            sol_rows = services.solvent.get_workbook_list(self.workbook)

        return sol_rows


def new():
    return json.dumps(
        {
            "amount_units": "mmol",
            "mass_units": "mg",
            "volume_units": "mL",
            "solvent_volume_units": "mL",
            "product_amount_units": "mmol",
            "product_mass_units": "mg",
            # reactants
            "reactant_masses": [],
            "reactant_names": [],
            "reactant_molecular_weights": [],
            "reactant_masses_raw": [],
            "reactant_amounts": [],
            "reactant_amounts_raw": [],
            "reactant_volumes": [],
            "reactant_volumes_raw": [],
            "reactant_equivalents": [],
            "reactant_physical_forms": [],
            "reactant_densities": [],
            "reactant_concentrations": [],
            "reactant_mns": [],
            # reagents
            "reagent_names": [],
            "reagent_molecular_weights": [],
            "reagent_densities": [],
            "reagent_concentrations": [],
            "reagent_amounts": [],
            "reagent_amounts_raw": [],
            "reagent_equivalents": [],
            "reagent_physical_forms": [],
            "reagent_hazards": [],
            "reagent_masses": [],
            "reagent_masses_raw": [],
            "reagent_volumes": [],
            "reagent_volumes_raw": [],
            # solvents
            "solvent_ids": [],
            "solvent_volumes": [],
            "solvent_names": [],
            "solvent_concentrations": [],
            "solvent_hazards": [],
            "solvent_physical_forms": [],
            # products
            "product_names": [],
            "product_mns": [],
            "product_amounts": [],
            "product_amounts_raw": [],
            "product_equivalents": [],
            "product_masses": [],
            "product_masses_raw": [],
            "product_physical_forms": [],
            # params
            "reaction_time": "",
            "reaction_atmosphere": "-select-",
            "reaction_temperature": "",
            "reaction_temperature_units": "",
            "reaction_time_units": "",
            "response_names": [],
            "response_units": [],
            "response_values": [],
            "response_notes": [],
        }
    )
