import re
from typing import List, Tuple

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

    def __init__(self, smiles, reaction_table_data, units):
        self.smiles = smiles
        self.previous_data = reaction_table_data
        self.units = units

    def update(self):
        pass

    def reload(self):
        pass

    def new(self):
        pass

    def render(self):
        pass
