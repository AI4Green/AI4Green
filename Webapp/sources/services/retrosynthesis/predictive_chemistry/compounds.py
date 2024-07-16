import itertools
import re
from typing import List, Tuple

import pandas as pd
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
        - Compound objects or "No Compound" if not found
        - Names of each compound or "No Compound" if not found
        - IDs of each compound or "No Compound" if not found
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
            compound_object = services.compound.get_compound_from_smiles(smiles)

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
        # same ion method as used elsewhere
        smiles_list = services.ions.update_ion_containing_list(
            compound_smiles.split(".")
        )
        # todo remove ion function if the services code works
        # older ion method written for this code
        # smiles_list = services.ions.update_ion_containing_list(
        #     compound_smiles.split(".")
        # )
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


class IonicCompounds:
    """
    Class to process ions in a smiles string. Tries to combine ions into an ionic compound. e.g., Na+ Cl- -> NaCl
    """

    def __init__(self, delimited_smiles: str):
        """
        Args:
            delimited_smiles: smiles with a '.' delimiter
        """

        self.smiles_string = delimited_smiles
        self.smiles_list = delimited_smiles.split(".")
        self.possible_ionic_compound_list = []
        self.ionic_compound_eval_list = []  # list of dicts
        self.ions_to_keep = []
        self.processed_smiles_list = []

    def get_smiles_list(self) -> List[str]:
        """
        Process ions in the smiles string.

        Returns:
            list: Processed smiles list.
        """
        self.make_ions_compounds()
        self.evaluate_ionic_compounds()
        self.select_best_ion_compounds()
        self.format_smiles_list()
        return self.processed_smiles_list

    def make_ions_compounds(self):
        """
        Make ionic compounds using itertools.
        """
        self.possible_ionic_compound_list = []
        for i in range(1, len(self.smiles_list) + 1):
            combinations = [
                y
                for y in [
                    x for x in itertools.combinations(enumerate(self.smiles_list), i)
                ]
            ]
            dic_ls = [{"idx": x[0], "ion": x[1]} for x in combinations[0]]
            ionic_compound = ""
            idx_ls = []
            for idx, dic in enumerate(dic_ls):
                ionic_compound += dic["ion"]
                idx_ls.append(dic["idx"])
            if "+" in ionic_compound or "-" in ionic_compound:
                new_dic = {"idx_list": idx_ls, "ionic_compound": ionic_compound}
                self.possible_ionic_compound_list.append(new_dic)

    def evaluate_ionic_compounds(self):
        """
        Evaluate ionic compounds. VH is likely to exist and L unlikely
        """
        for ionic_compound_dict in self.possible_ionic_compound_list:
            ion_smiles = ionic_compound_dict["ionic_compound"]
            compound_object = services.compound.get_compound_from_smiles(ion_smiles)
            balanced_charge = self.charges_balanced(ion_smiles)
            if compound_object and balanced_charge:
                ionic_compound_dict.update({"eval": "H"})
            elif compound_object:
                ionic_compound_dict.update({"eval": "M"})
            elif balanced_charge:
                ionic_compound_dict.update({"eval": "L"})
            else:
                ionic_compound_dict.update({"eval": "VL"})
            self.ionic_compound_eval_list.append(ionic_compound_dict)

    def select_best_ion_compounds(self):
        """
        Select the best ion compounds.
        """
        df = pd.DataFrame(self.ionic_compound_eval_list)
        score_list = ["H", "M", "L", "VL"]
        for score in score_list:
            best_ion = df[df["eval"] == score]
            if not best_ion.empty:
                self.ions_to_keep.append(best_ion.to_dict(orient="records")[0])
                break

    def format_smiles_list(self):
        """
        Format the smiles list.
        """
        ions_in_previous_salts = 0
        self.processed_smiles_list = self.smiles_list
        for ion_dict in self.ions_to_keep:
            ion_indexes = ion_dict["idx_list"]
            insert_index = ion_indexes[0] - ions_in_previous_salts
            end_index = insert_index + len(ion_indexes)
            del self.processed_smiles_list[insert_index:end_index]
            self.processed_smiles_list.insert(insert_index, ion_dict["ionic_compound"])

    @staticmethod
    def charges_balanced(ion_string: str) -> bool:
        """
        Check if charges are balanced in the ion string.

        Returns:
            bool: True if charges are balanced, False otherwise.
        """
        positive_charges = ion_string.count("+")
        negative_charges = ion_string.count("-")
        if positive_charges == negative_charges:
            return True
        return False


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
