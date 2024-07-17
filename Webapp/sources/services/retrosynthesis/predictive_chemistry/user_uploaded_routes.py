import base64
import io
import uuid
from typing import Dict, List, Tuple, Union

import numpy as np
import pandas as pd
from sources import services

from .compounds import get_compound_data, ion_check
from .utils import sig_figs_on_numbers


def read_user_route_file(contents: str, filename: str) -> Union[Tuple[Dict, Dict], str]:
    """
    Reads in a user supplied route file rather than a predicted retrosynthesis
    Args:
        contents - the file contents with the route details
        filename - the name of the file the user uploaded
    Returns:
        a dictionary with route data
        a dictionary with conditions data
    """
    content_type, content_string = contents.split(",")
    decoded = base64.b64decode(content_string)

    if "csv" in filename:
        route_df = pd.read_csv(io.StringIO(decoded.decode("utf-8")))
    elif "xls" in filename or "ods" in filename:
        route_df = pd.read_excel(io.BytesIO(decoded))
    else:
        return "error processing file. Did you upload a csv, xls, or ods file?"
    route_df = process_optional_columns(route_df)
    # make the routes dictionary
    route_object = RetroRoute(route_df)
    unique_identifier = str(uuid.uuid4())
    route_data, condition_data = route_object.to_dict_arrays()
    route_dict = {
        "routes": {"Route 1": {"steps": route_data, "score": 1}},
        "uuid": unique_identifier,
    }
    # make the conditions dictionary
    conditions_dict = {
        "routes": {"Route 1": condition_data},
        "uuid": unique_identifier,
    }
    return route_dict, conditions_dict


def process_optional_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    If these are present keep the values, if not present add default values
    Args:
        df - the dataframe under construction with the route/conditions data
    Returns:
         the updated dataframe
    """
    optional_string_columns = ["reagents", "catalyst", "solvents"]
    df[optional_string_columns] = df[optional_string_columns].replace(np.nan, "")
    return df


class Conditions:
    """A class to make conditions for a user uploaded route"""

    def __init__(
        self, step_dataframe: pd.DataFrame, reactants_smiles: str, product_smiles: str
    ):
        """Creates a new instance of the conditions class"""
        self.step_df = step_dataframe
        self.predicted_conditions = {}
        self.reactants_smiles = reactants_smiles
        self.product_smiles = product_smiles

    def to_dicts(self) -> Dict:
        step_df = self.step_df
        conditions_set = {
            "temperature": step_df[step_df["role"] == "temperature"]["value"].squeeze(),
            "solvent": step_df[step_df["role"] == "solvents"]["value"]
            .squeeze()
            .replace(";", "."),
            "catalyst": step_df[step_df["role"] == "catalyst"]["value"]
            .squeeze()
            .replace(";", "."),
            "reagents": step_df[step_df["role"] == "reagents"]["value"]
            .squeeze()
            .replace(";", "."),
        }
        self.predicted_conditions = {"Condition Set 1": conditions_set}
        self.process_conditions()
        return self.predicted_conditions

    def to_dicts2(self):
        step_df = self.step_df
        conditions_set = {
            "temperature": step_df["temperature"],
            "solvent": step_df["solvents"].replace(";", "."),
            "catalyst": step_df["catalyst"].replace(";", "."),
            "reagents": step_df["reagents"].replace(";", "."),
        }
        self.predicted_conditions = {"Condition Set 1": conditions_set}
        self.process_conditions()
        return self.predicted_conditions

    def process_conditions(self):
        """Constructs a list of dictionaries with keys for the Names and Database IDS of:
        solvents, reagent, catalyst, reactants, and products, with the properties as a list
        """
        compound_keys = ["solvent", "reagents", "catalyst", "reactant", "product"]
        for condition_set in self.predicted_conditions.values():
            # add all keys to the condition_set dict, so it has name/id keys for all compound types.
            self.add_compound_to_dict(self.reactants_smiles, condition_set, "reactant")
            self.add_compound_to_dict(self.product_smiles, condition_set, "product")
            self.add_name_keys_to_dict(condition_set, compound_keys)
            for compound_type, compound_smiles in condition_set.items():
                # skips non-compound keys in condition set
                if compound_type not in compound_keys or not isinstance(
                    compound_smiles, str
                ):
                    continue
                # get list of smiles then compounds, then save name and ids to the dictionary
                smiles_list, compound_list, name_list, id_list = get_compound_data(
                    compound_smiles, compound_type
                )
                self.update_conditions_dict(
                    name_list, id_list, smiles_list, condition_set, compound_type
                )
                sig_figs_on_numbers(condition_set)

    @staticmethod
    def add_compound_to_dict(smiles, condition_set, compound_type):
        """Adds smiles as a string delimited by '.'"""
        if type(smiles) is str:
            condition_set[compound_type] = smiles
        elif type(smiles) is list:
            try:
                condition_set[compound_type] = ".".join(smiles)
            except Exception as e:
                print(e)
        else:
            print("unexpected type")

    @staticmethod
    def add_name_keys_to_dict(condition_set, compound_keys):
        for key in compound_keys:
            condition_set[f"{key}_names"] = ""
            condition_set[f"{key}_ids"] = ""
            condition_set[f"{key}_smiles"] = ""

    @staticmethod
    def update_conditions_dict(
        name_list, id_list, smiles_list, condition_set, compound_type
    ):
        # add name and id for each compound
        condition_set[f"{compound_type}_names"] = name_list
        condition_set[f"{compound_type}_ids"] = id_list
        condition_set[f"{compound_type}_smiles"] = smiles_list

    def smiles_str_to_list(self, compound_smiles: str):
        """Takes the smiles string with '.' delimiter and makes a list. Checking for ions and handling them"""
        if ion_check(compound_smiles):
            # function to process ions then convert to list, eg na+.oh-.cco > [na.oh, cco]
            smiles_list = services.retrosynthesis.predictive_chemistry.compounds.smiles_str_to_list(
                compound_smiles
            )
            print(smiles_list)
        else:
            print(compound_smiles)
            smiles_list = compound_smiles.split(".")
        return smiles_list

    def smiles_list_to_compounds(self, smiles_list: List):
        compound_ls = []
        for smiles in smiles_list:
            compound = services.compound.get_compound_from_smiles(smiles)
            compound_ls.append(compound)
        return compound_ls

    @staticmethod
    def get_compound_name_list(compound_list, compound_type):
        name_list = []
        for compound in compound_list:
            if compound == "Not Found":
                name = "Not Found"
            elif compound == "No Compound":
                name = "No Compound"
            else:
                if compound_type == "solvent":
                    try:
                        name = compound.solvent.name
                    except AttributeError:
                        name = compound.name
                        print(name)
                else:
                    name = compound.name
            name_list.append(name)
        return name_list

    @staticmethod
    def get_compound_id_list(compound_list, compound_type):
        id_list = []
        for compound in compound_list:
            if compound == "Not Found":
                compound_id = "Not Found"
            elif compound == "No Compound":
                compound_id = "No Compound"
            else:
                compound_id = compound.id
            id_list.append(compound_id)
        return id_list


class RetroRoute:
    """For a route dictionary, extracts the reaction routes"""

    def __init__(self, route_dataframe):
        self.reactions = []
        self.depth = 0
        self.previous_depth = 0
        self.route_df = route_dataframe
        self.parent = None
        self.parent_smiles = None
        self.node_id_ls = []
        self.retrosynthesis_nodes = []
        self.condition_sets = {}

    def to_dict_arrays(self):
        route_df = self.route_df
        products = route_df["product"]
        # go through each product and process each step into the data format used
        for product in products:
            step_df = route_df[route_df["product"] == product].squeeze()
            reactants = step_df["reactants"].split(";")
            # reactants_list = reactants.split(';')

            # retrosynthesis
            # self.make_node_dict(step_df)
            self.make_node_dict2(step_df)
            # conditions
            conditions_object = Conditions(step_df, reactants, product)
            self.condition_sets.update(
                {self.node_id_ls[-1]: conditions_object.to_dicts2()}
            )
            print(type(reactants), reactants)
            terminal_node_reactants = list(set(reactants) - set(products))
            if terminal_node_reactants:
                for reactant in terminal_node_reactants:
                    self.make_terminal_node_dict(reactant)

            # condition_step = format_condition_data(step_df)
        print(self.node_id_ls)
        return self.retrosynthesis_nodes, self.condition_sets

        # self.add_terminal_nodes()

    def make_terminal_node_dict(self, reactant):
        node_dict = {
            "node_id": self.node_id(),
            "smiles": reactant,
            "node": None,
            "depth": self.depth,
            "parent": self.parent,
            "parent_smiles": self.parent_smiles,
            "child_smiles": [],
            "reaction_class": [],
        }
        self.retrosynthesis_nodes.append(node_dict)

    def make_node_dict(self, step_df):
        """ """

        node_dict = {
            "node_id": self.node_id(),
            "smiles": step_df[step_df["role"] == "product"]["value"].squeeze(),
            "node": None,
            "depth": self.depth,
            "parent": self.parent,
            "parent_smiles": self.parent_smiles,
            "child_smiles": step_df[step_df["role"] == "reactants"]["value"]
            .squeeze()
            .split(";"),
            "reaction_class": step_df.get("reaction_class"),
        }
        print(node_dict)
        self.retrosynthesis_nodes.append(node_dict)
        self.parent = node_dict
        self.parent_smiles = node_dict["smiles"]
        self.depth += 1

    def make_node_dict2(self, step_df):
        """ """
        print()

        node_dict = {
            "node_id": self.node_id(),
            "smiles": step_df["product"],
            "node": None,
            "depth": self.depth,
            "parent": self.parent,
            "parent_smiles": self.parent_smiles,
            "child_smiles": step_df["reactants"].split(";"),
            "reaction_class": step_df.get("reaction_class"),
        }
        print(node_dict)
        self.retrosynthesis_nodes.append(node_dict)
        self.parent = node_dict
        self.parent_smiles = node_dict["smiles"]
        self.depth += 1

    def node_id(self):
        """Creates an id for the node"""
        if self.node_id_ls:
            id_ls = [x.split("node-")[-1] for x in self.node_id_ls]
            # 0-1
            current_depth_id_ls = [
                int(x.split(f"{self.depth}-")[-1])
                for x in id_ls
                if f"{self.depth}-" in x
            ]
            if current_depth_id_ls:
                next_number = max(current_depth_id_ls) + 1
            else:
                next_number = 0
            # e.g., node-0, node-1-0
            node_id = f"node-{self.depth}-{next_number}"
        else:
            # for first node
            node_id = "node-0"
        self.node_id_ls.append(node_id)
        return node_id
