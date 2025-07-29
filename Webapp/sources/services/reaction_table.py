import re
from typing import List, Tuple

from flask import json, jsonify, render_template
from sources import services


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

        # reaction table species
        self.reactants = []
        self.products = []
        self.reagents = []
        self.solvents = []

        # interactive elements
        self.solvent_dropdown = self.get_solvent_dropdown()

        # load reaction_table from provided reaction
        self.load()

    def load(self):
        compounds, units = services.compound.SketcherCompound.from_reaction_table_dict(
            json.loads(self.reaction.reaction_table_data), self.workbook
        )
        self.units = units
        self.reactants = compounds["reactant"]
        self.products = compounds["product"]
        self.reagents = compounds["reagent"]
        self.solvents = compounds["solvent"]

    def update(self, smiles):
        pass

    def reload(self):
        pass

    def new(self):
        pass

    def render(self):
        reaction_table = render_template(
            "reactions/_reaction_table.html",
            reactants=self.reactants,
            reagents=self.reagents,
            solvents=self.solvents,
            products=self.products,
            number_of_reactants=len(self.reactants),
            number_of_reagents=len(self.reagents),
            number_of_solvents=len(self.solvents),
            number_of_products=len(self.products),
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
        }

    def get_solvent_dropdown(self):
        if self.demo == "demo":
            sol_rows = services.solvent.get_default_list()
        else:
            sol_rows = services.solvent.get_workbook_list(self.workbook)

        return sol_rows
