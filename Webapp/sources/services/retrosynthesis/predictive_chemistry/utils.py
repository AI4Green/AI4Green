import base64
from typing import Dict, List, Literal, Union

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

from ..style_sheets import CytoscapeStyles as cytoStyles

scale_factor = cytoStyles.scale_factor


def sig_figs_on_numbers(condition_set: Dict):
    """
    Round the values for score and temperature to the specified number of significant figures.
    Args:
        condition_set - the condition set we are rounding values for.

    """
    keys_to_replace = {"score": 2, "temperature": 0}
    for key, value in condition_set.items():
        if key in keys_to_replace.keys():
            condition_set[key] = round(value, keys_to_replace[key])


def rdkit_smiles_to_image(smiles: str):
    """Single molecule image
    Args:
        smiles - SMILES string for the compound we are making an image of.
    """
    mol = Chem.MolFromSmiles(smiles)
    d = rdMolDraw2D.MolDraw2DCairo(350, 300)
    d.DrawMolecule(mol)
    d.FinishDrawing()
    img_data = d.GetDrawingText()
    img_data = base64.b64encode(img_data)
    img_data = img_data.decode()
    img_data = "{}{}".format("data:image/png;base64,", img_data)
    return img_data


def reaction_smiles_to_image(window_width: int, smiles: str) -> str:
    """
    Take a reaction smiles string and returns an image data string
    Args:
        window_width - the width of the user's window used to scale images to an appropriate size
        smiles - the reaciton smiles we are making an image of.
    Returns:
        base64 png image string to show in HTML
    """
    rxn = AllChem.ReactionFromSmarts(smiles, useSmiles=True)
    width = round(window_width / 3.4)
    height = round(width / 2.5)
    drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
    drawer.drawOptions().padding = 0.0
    drawer.SetFontSize(6)
    drawer.maxFontSize = 6
    drawer.DrawReaction(rxn)
    drawer.FinishDrawing()
    image_file = drawer.GetDrawingText()
    img_data = base64.b64encode(image_file)
    img_data = img_data.decode()
    img_data = "{}{}".format("data:image/png;base64,", img_data)
    return img_data


def functionality_disabled_check(
    functionality_status: Union[Literal["disabled"], Literal["enabled"]]
) -> bool:
    """
    Many functionalities are disabled if a user is not in a workbook.
    Args:
        functionality_status - either disabled or enabled

    Returns - True if the functionality is disabled.
    """
    if functionality_status == "disabled":
        return True
    return False


def get_current_route(routes: List[Dict], selected_route: str) -> Dict:
    """
    Gets the current route from the list of routes
    Args:
        routes - the list of routes
        selected_route - the selected route in the form e.g., 'Route 1'
    Returns:
        The selected route as a dictionary
    """
    selected_route_idx = int(selected_route[-1]) - 1  # minus 1 to use zero-based index
    current_route = routes[selected_route_idx]
    return current_route


def encodings_to_smiles_symbols(input_str: str) -> str:
    """
    Returns the string with the SMILES symbols restored, replacing the URL encodings.
    Args:
        input_str - the SMILES string being decoded
    Returns:
        the decoded smiles string

    """
    return (
        input_str.replace("%23", "#")
        .replace("%2B", "+")
        .replace("%2D", "-")
        .replace("%40", "@")
    )


def smiles_symbols_to_encodings(input_str: str) -> str:
    """
    Returns the string with the URL encodings, replacing the SMILES symbols.
    Args:
        input_str - the SMILES string being encoded
    Returns:
        the encoded SMILES string
    """
    return (
        input_str.replace("#", "%23")
        .replace("+", "%2B")
        .replace("-", "%2D")
        .replace("@", "%40")
    )


def smiles_not_valid(smiles_regex: str) -> bool:
    """
    Returns a bool to indicate if SMILEs is valid.
    Args:
        smiles_regex - the regex used to detect if SMILES is valid, from the field in the frontend
    Returns:
        True if the smiles is invalid
    """
    if "_invalid" in smiles_regex:
        return True
    return False
