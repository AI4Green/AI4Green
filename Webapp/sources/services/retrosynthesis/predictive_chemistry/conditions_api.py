import json
from typing import Dict, List, Tuple, Union

import requests
from flask import current_app

from .compounds import get_compound_data
from .utils import encodings_to_smiles_symbols, sig_figs_on_numbers


def get_conditions(solved_routes: dict) -> Dict[str, dict]:
    """
    Takes the retrosynthetic routes and send an api call with the reaction smiles
    to get conditions for each forward reaction.

    Args:
        solved_routes - retrosynthetic routes

    Returns:
        List of conditions with a dictionary for each route.

    """
    conditions_results = {}
    for (
        route_label,
        route,
    ) in solved_routes.items():
        conditions_results.update(
            RouteConditions(route_label, route).get_route_conditions()
        )
    return conditions_results


class RouteConditions:
    """
    Class to handle conditions for a specific route.
    """

    base_url = current_app.config["CONDITIONS_API_URL"]
    api_key = current_app.config["CONDITIONS_API_KEY"]

    def __init__(self, route_label: str, route: dict):
        self.route_label = route_label
        self.route = route

    def get_route_conditions(self) -> Dict:
        """
        Get conditions for each step in the route.

        Returns:
            dict: Conditions for the route.
        """
        route_conditions = {}
        for node in self.route["steps"]:
            # child smiles are the reactants. if there are none it is a terminal node.
            if self._not_terminal(node):
                # get the conditions or make note of the failed api call
                api_status, reaction_conditions = ReactionConditions(
                    node, self
                ).get_reaction_conditions()

                if api_status == "failed":
                    reaction_conditions = "Condition Prediction Unsuccessful"

                route_conditions.update({node["node_id"]: reaction_conditions})

        return {self.route_label: route_conditions}

    @staticmethod
    def _not_terminal(node) -> bool:
        """Returns True if the node has child smiles"""
        return node.get("child_smiles")


class ReactionConditions(RouteConditions):
    """
    Class to handle conditions for a specific reaction.
    """

    def __init__(self, node, route_conditions):
        super().__init__(route_conditions.base_url, route_conditions.api_key)
        self.reactants_smiles = node["child_smiles"]
        self.product_smiles = node["smiles"]
        self.reaction_smiles = encodings_to_smiles_symbols(
            ".".join(self.reactants_smiles) + ">>" + self.product_smiles
        )
        self.predicted_conditions = None
        self.api_status = ""

    def get_reaction_conditions(self) -> Tuple[str, Dict]:
        """
        Gets conditions for a specific reaction

        Returns
            dict - the api status success or fail and why if failure
            list of dicts - predicted conditions

        """
        self._api_call()
        processed_conditions = ProcessConditions(
            self.predicted_conditions, self.reactants_smiles, self.product_smiles
        ).process_conditions()
        return self.api_status, processed_conditions

    def _api_call(self):
        """
        Makes the api call to the condition prediction url to get the conditions data
        """
        url = f"{self.base_url}/api/v1/condition_cleaned/"
        smiles = self.reaction_smiles
        data = {"smiles": smiles, "n_conditions": 10}

        # Make the POST request
        try:
            response = requests.post(
                url, json=data, headers={"Content-Type": "application/json"}
            )
            self.predicted_conditions = self._process_api_response(response)
            self.api_status = "success"
        except requests.exceptions.HTTPError as e:
            print("Request failed", e, url)
            self.api_status = "failed"
            self._exit(e, "")

    @staticmethod
    def _process_api_response(response) -> dict:
        """Labels the condition sets from the API response"""
        conditions = json.loads(response.content)
        conditions_dict = {}
        for idx, condition_set in enumerate(conditions, 1):
            conditions_dict.update(
                {
                    f"Condition Set {idx}": {
                        "temperature": condition_set[0],
                        "solvent": condition_set[1],
                        "reagents": condition_set[2],
                        "catalyst": condition_set[3],
                    }
                }
            )
        return conditions_dict

    def _exit(self, error: requests.exceptions.HTTPError, conditions: str):
        return self.api_status, error, conditions


class ProcessConditions:
    """
    A Class to process conditions got from the api request into a list.
    Processing includes looking up compounds in the database to get additional data and rounding numbers.
    """

    def __init__(
        self, predicted_conditions: dict, reactants_smiles: str, product_smiles: str
    ):
        self.predicted_conditions = predicted_conditions
        self.reactants_smiles = reactants_smiles
        self.product_smiles = product_smiles

    def process_conditions(self) -> Dict[str, Dict]:
        """Constructs a dictionary of dictionaries with keys for the Names and Database IDS of:
        solvents, reagent, catalyst, reactants, and products, with the properties as a list
        """
        compound_keys = ["solvent", "reagents", "catalyst", "reactant", "product"]
        for condition_set in self.predicted_conditions.values():
            # add all keys to the condition_set dict, so it has name/id keys for all compound types.
            self._add_compound_to_dict(self.reactants_smiles, condition_set, "reactant")
            self._add_compound_to_dict(self.product_smiles, condition_set, "product")
            self._add_name_keys_to_dict(condition_set, compound_keys)
            for compound_type, compound_smiles in condition_set.items():
                # skips non-compound keys in condition set
                if compound_type not in compound_keys:
                    continue

                # get list of smiles then compounds, then save name and ids to the dictionary
                smiles_list, compound_list, name_list, id_list = get_compound_data(
                    compound_smiles, compound_type
                )

                self._update_conditions_dict(
                    name_list, id_list, smiles_list, condition_set, compound_type
                )
                sig_figs_on_numbers(condition_set)
        return self.predicted_conditions

    @staticmethod
    def _add_compound_to_dict(
        smiles: Union[str, List], condition_set: Dict, compound_type: str
    ):
        """
        Adds smiles as a string delimited by '.' or list
        Args:
            smiles - the SMILES string being added as either a delimited string or a list
            condition_set - the dictionary we are adding the SMILES to.
            compound_type - whether the SMILES are reactants, reagents, solvents, catalysts, or products.
        """
        if type(smiles) is str:
            condition_set[compound_type] = smiles
        elif type(smiles) is list:
            try:
                condition_set[compound_type] = ".".join(smiles)
            except Exception as e:
                print(e)
        else:
            raise TypeError("SMILES to add to condition set is not type list or string")

    @staticmethod
    def _add_name_keys_to_dict(condition_set: dict, compound_keys: List[str]):
        """
        Initialises keys for the condition set dictionary
        Args:
            condition_set - the dictionary we are initialising keys for.
            compound_keys - whether the compound is type reactant, reagent, solvent, product.
        """
        for key in compound_keys:
            condition_set[f"{key}_names"] = ""
            condition_set[f"{key}_ids"] = ""
            condition_set[f"{key}_smiles"] = ""

    @staticmethod
    def _update_conditions_dict(
        name_list: List[str],
        id_list: List[int],
        smiles_list: List[str],
        condition_set: Dict,
        compound_type: str,
    ):
        """
        Adds the names, ids, and SMILEs for each compound to the dictionary
        Args:
            name_list - the list of compound names
            id_list - the list of primary key database ids of the compounds
            smiles_list - the list of SMILES strings of the compounds
            condition_set - the dictionary being updated
            compound_type - whether the compounds are reactants, reagents, solvents, products, or catalysts.
        """
        condition_set[f"{compound_type}_names"] = name_list
        condition_set[f"{compound_type}_ids"] = id_list
        condition_set[f"{compound_type}_smiles"] = smiles_list
