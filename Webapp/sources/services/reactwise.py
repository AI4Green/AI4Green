from typing import Dict, List, Union

import requests
from flask import current_app


class ReactWiseStep:

    """Class that processes and contains info from reactwise steps"""

    def __init__(self, step_id: int):
        self.step_id = step_id
        self.categorical_inputs = None
        self.categorical_map = {}

    def load_step_experiments(self):
        inputs_url = f"https://api.reactwise.com/views/step-experiments/{self.step_id}"  # noqa: E231
        components_url = (
            f"https://api.reactwise.com/step-components/{self.step_id}"  # noqa: E231
        )
        print(components_url)

        headers = {
            "Authorization": "Bearer " + current_app.config["REACTWISE_API_KEY"]
        }  # possibly add reactwise key per user?

        input_response = requests.request(
            "GET", inputs_url, headers=headers
        ).json()  # not working yet

        # component_response = requests.request(
        #     "GET", components_url, headers=headers
        # ).json()

        self._extract_categorical(input_response)
        self._extract_continuous(input_response)
        self._extract_experiment_inputs(input_response)

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

    print(sorted_steps)
    return sorted_steps

    # for step in step_response.get("data", []):
    #     print(step)
