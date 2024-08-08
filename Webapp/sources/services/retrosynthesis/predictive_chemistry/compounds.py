import re
from typing import List, Tuple

from sources import services

"""
File containing functions and classes that get data on compounds
and interconvert between different molecular representations
"""


def get_compound_data(compound_smiles: str, role: str) -> Tuple[List, List, List, List]:
    """
    Takes the list of smiles and type of compound to look up the compounds in the database and return data on them.

    Args:
        compound_smiles - '.' delimited string of the SMILES of each compound of a particular type (e.g., reactants)
        role - the type of compound ina  reaction - reactant, solvent, etc.

    Returns:
        A tuple containing lists of:
        - SMILES for each compound
        - Compound objects or "Not Found" if not found
        - Names of each compound or "Not Found" if not found
        - IDs of each compound or "Not Found" if not found
    """
    smiles_list = smiles_str_to_list(compound_smiles)
    compound_data = {"smiles": [], "compound_object": [], "name": [], "id": []}
    for smiles in smiles_list:
        # if there are none of a role in a reaction, no solvents for example
        if not smiles:
            compound_data["compound_object"].append(None)
            compound_data["name"].append(f"No {role}")
            compound_data["id"].append(None)
            compound_data["smiles"].append(None)
        else:
            compound_object = services.compound.from_smiles(smiles)

            compound_data["smiles"].append(smiles)
            if compound_object:
                compound_data["compound_object"].append(compound_object)
                compound_data["name"].append(compound_object.name)
                compound_data["id"].append(compound_object.id)
            else:
                # if there is a compound, but it is not in the AI4Green database
                compound_data["compound_object"].append("Not Found")
                compound_data["name"].append("Not Found")
                compound_data["id"].append("Not Found")

    return (
        compound_data["smiles"],
        compound_data["compound_object"],
        compound_data["name"],
        compound_data["id"],
    )


def smiles_str_to_list(compound_smiles: str) -> List[str]:
    """
    Takes the smiles string with '.' delimiter and makes a list. Checking for ions and handling them

    Args:
        compound_smiles (str): The SMILES string with '.' delimiter.

    Returns:
        A list of SMILES strings, processed to handle ions if present.
    """

    if ion_check(compound_smiles):
        if "O=C([O-])[O-]" in compound_smiles:
            compound_smiles = fix_carbonates(compound_smiles)
        smiles_list = services.ions.update_ion_containing_list(
            compound_smiles.split(".")
        )
    else:
        smiles_list = compound_smiles.split(".")
    return smiles_list


def fix_carbonates(compound_smiles: str) -> str:
    """
    Fixes an issue where carbonates are returned from the conditions missing one of the counterions
    For example, NaCO3  --> Na2CO3 and KCO3 -> K2CO3 for sodium/potassium carbonate
    Args:
        compound_smiles - the SMILES string for the compound being fixed
    Returns:
        The reformed SMILES string to return a neutral carbonate species
    """
    carbonate = "O=C([O-])[O-]"
    compounds = compound_smiles.split(".")

    previous_compound = ""
    for idx, compound in enumerate(compounds):
        if previous_compound == carbonate:
            counter_ion = compound
            # Add a second counter ion to the string after the ion pattern
            compounds.insert(idx + 1, counter_ion)
            break
        previous_compound = compound

    # Join the modified compounds back into a single SMILES string
    new_compound_smiles = ".".join(compounds)
    return new_compound_smiles


def ion_check(smiles: str) -> bool:
    """
    Checks whether a compound contains ions.
    Args:
        smiles - the SMILES string being checked.

    Returns:
        True if the compound contains ions.
    """
    if re.findall(r"(\+]|-]|\b\d+\])", smiles):
        return True
    return False
