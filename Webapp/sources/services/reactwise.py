import json
from typing import Dict, List, Union

import requests
from flask import current_app
from sources import services


class ReactWiseStep:

    """Class that processes and contains info from reactwise steps"""

    def __init__(self, step_id: int):
        self.step_id = step_id
        self.categorical_inputs = None
        self.continuous_map = None
        self.categorical_map = None
        self.experimental_details = None
        self.reaction_smiles = None
        self.reactants = None
        self.products = None

        self.load_step_experiments()

    def _get_api_call(self, url):
        headers = {
            "Authorization": "Bearer " + current_app.config["REACTWISE_API_KEY"],
        }  # possibly add reactwise key per user?

        payload = json.dumps(
            {
                "step_id": str(self.step_id),
            }
        )
        # currently only get requests
        return requests.request("GET", url, headers=headers, data=payload).json()

    def load_step_experiments(self):
        input_response = self._get_api_call(
            f"https://api.reactwise.com/views/step-experiments/{self.step_id}"  # noqa: E231
        )
        component_response = self._get_api_call(
            "https://api.reactwise.com/step-components"  # noqa: E231
        )

        self._extract_categorical(input_response)
        self._extract_continuous(input_response)
        self._extract_experiment_inputs(input_response)
        self._extract_reaction_smiles(component_response)

    def _extract_categorical(self, response: Dict[str, Union[int, str]]):
        """
        Converts reactwise categorical inputs to simple lookup dict
        """
        categorical_inputs = response.get("categorical_inputs", {})
        categorical_input_values = response.get("categorical_input_values", {})

        # dictionary comprehension to handle nested dicts
        self.categorical_map = {
            input_item.get("id"): {
                "name": input_item.get("name"),
                "values": {
                    val.get("id"): val.get("value")
                    for val in categorical_input_values.get(
                        str(input_item.get("id")), []
                    )
                },
            }
            for input_item in categorical_inputs
        }

    def _extract_continuous(self, response: Dict[str, Union[int, str]]):
        continuous_inputs = response.get("continuous_inputs", {})
        self.continuous_map = {
            input_item.get("id"): {
                "name": input_item.get("name"),
                "unit": input_item.get("unit"),
            }
            for input_item in continuous_inputs
        }

    def _extract_experiment_inputs(self, response: Dict[str, Union[int, str]]):
        experiment_details = {}

        categorical_inputs = response.get("experiment_categorical_inputs", {})
        continuous_inputs = response.get("experiment_continuous_inputs", {})

        for exp_id, data in categorical_inputs.items():
            exp_dict = {}
            exp_continuous = continuous_inputs.get(exp_id)
            for inp in data:
                input_id = inp.get("categorical_input_id")
                value_id = inp.get("categorical_input_value_id")

                map_entry = self.categorical_map.get(input_id)

                name = map_entry.get("name")
                value = map_entry.get("values").get(value_id)

                exp_dict[name] = value

            for cont_inp in exp_continuous:
                exp_dict[
                    self.continuous_map.get(cont_inp.get("continuous_input_id")).get(
                        "name"
                    )
                ] = {
                    "value": cont_inp.get("value"),
                    "unit": self.continuous_map.get(
                        cont_inp.get("continuous_input_id")
                    ).get("unit"),
                }
            experiment_details[exp_id] = exp_dict
        self.experimental_details = experiment_details

    def _extract_reaction_smiles(self, api_response):
        """
        Possibly need to handle multiple reaction smiles for different reactions in a step?
        """
        step_components = [
            x for x in api_response["data"] if x.get("step_id", None) == self.step_id
        ]

        reactant_smiles = []
        product_smiles = []

        for component in step_components:
            step_component = self._get_api_call(
                f"https://api.reactwise.com/components/{component.get('component_id')}"  # noqa: E231
            )
            component_smiles = step_component.get("data").get("smiles")
            if component.get("reactant"):
                reactant_smiles.append(component_smiles)
            else:
                product_smiles.append(component_smiles)

        self.reactants = reactant_smiles
        self.products = product_smiles
        self.reaction_smiles = (
            f"{'.'.join(reactant_smiles)}>>{'.'.join(product_smiles)}"
        )


def list_steps():
    """
    Lists all reactwise steps arranged by project
    """
    projects_url = "https://api.reactwise.com/projects"
    steps_url = "https://api.reactwise.com/steps"

    headers = {
        "Authorization": "Bearer " + current_app.config["REACTWISE_API_KEY"]
    }  # possibly add reactwise key per user?

    projects = (
        requests.request("GET", projects_url, headers=headers, data={})
        .json()
        .get("data", None)
    )

    steps = (
        requests.request("GET", steps_url, headers=headers, data={})
        .json()
        .get("data", None)
    )

    sorted_steps = {}

    for project in projects:
        pid = project.get("id", None)
        psteps = [x for x in steps if x.get("project_id") == pid]

        sorted_steps[project.get("name")] = psteps
    return sorted_steps


# should write an individual function to extract each recognised element and add to extract function
def assign_reaction_time(reactwise_experiment_details, reaction_table_data):
    time_info = reactwise_experiment_details.get("Time", {})
    (
        reaction_table_data["reaction_time"],
        reaction_table_data["reaction_time_units"],
    ) = time_info.get("value"), time_info.get("unit")


def assign_reaction_temperature(reactwise_experiment_details, reaction_table_data):
    temp_info = reactwise_experiment_details.get("Temperature", {})
    (
        reaction_table_data["reaction_temperature"],
        reaction_table_data["reaction_temperature_units"],
    ) = temp_info.get("value"), temp_info.get("unit")


def extract_reaction_solvents(reactwise_experiment_details, unknown_solvents):
    """
    Gets reaction solvent by name and returns the solvent db object. If name is not in db, only the name is returned
    """
    solvent_name = reactwise_experiment_details.get("Solvent", "")
    solvent = services.solvent.from_name(solvent_name)

    if solvent:
        return solvent
    else:
        unknown_solvents.append(solvent_name)
        return None


def extract_temp_time_fields(reactwise_experiment_details, reaction_table_data):
    assign_reaction_time(reactwise_experiment_details, reaction_table_data)
    assign_reaction_temperature(reactwise_experiment_details, reaction_table_data)
