import pprint
from typing import Dict, List, Union

import requests
from flask import current_app, render_template, request
from flask_login import current_user
from sources import models, services

from . import reaction_set_bp


@reaction_set_bp.route("/reaction_set")
def reaction_set():
    # to do
    # add workbook, workgroup to reaction set page
    # fix atuosave sketcher to only update reaction table in set mode (do we want to try autosaving?)
    # then fix apply to well and apply to all
    # colours for unsaved/edited wells
    # import from reactwise (try at least one small bit of data)
    #
    return render_template("reaction_set.html")


@reaction_set_bp.route("/click_and_drag")
def click_and_drag():
    return render_template("click-and-drag.html")


@reaction_set_bp.route("/import_from_reactwise", methods=["GET", "POST"])
def import_from_reactwise():
    step_id = request.json.get("reactwiseID", None)
    workbook_name = request.json.get("workbook", None)
    workgroup_name = request.json.get("workgroup", None)

    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )

    data = ReactWiseStep(step_id)
    data.load_step_experiments()
    creator = services.user.person_from_current_user()

    set_name = "ReactWise Set " + str(step_id)
    reactions = []

    for reactwise_id, details in data.experimental_details.items():
        reaction_id = services.reaction.get_next_reaction_id_for_workbook(workbook.id)
        print(reaction_id)
        reaction = services.reaction.add(
            name="reactwise-" + reactwise_id,
            creator=creator,
            reaction_id=reaction_id,
            workbook_id=workbook.id,
            reaction_table={},
            summary_table={},
            reaction_smiles="CC.CC>>CC.CCC",  # change this obvs
        )
        reactions.append(reaction)

    set_id = services.reaction_set.next_id_in_workbook(workbook.id)
    reaction_set = services.reaction_set.add(
        name=set_name,
        set_id=set_id,
        creator=creator,
        workbook=workbook,
        reactions=reactions,
    )

    return render_template(
        "reaction_set.html",
        reaction_set=reaction_set,
    )


class ReactWiseStep:
    """Class that processes and contains info from reactwise steps"""

    def __init__(self, step_id: int):
        self.step_id = step_id
        self.categorical_inputs = None
        self.categorical_map = {}

    def load_step_experiments(self):
        inputs_url = f"https://api.reactwise.com/views/step-experiments/{self.step_id}"
        # components_url = f"https://api.reactwise.com/step-components/{self.step_id}"

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
