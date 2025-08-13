from typing import Dict, List, Union

import requests
from bokeh.embed import components
from flask import (
    Response,
    current_app,
    jsonify,
    redirect,
    render_template,
    request,
    url_for,
)
from flask_login import current_user, login_required
from sources import models, services

from . import reaction_set_bp


@reaction_set_bp.route("/reaction_set/<workgroup_name>/<workbook_name>/<set_id>")
@login_required
def reaction_set(workgroup_name, workbook_name, set_id):
    r_set = services.reaction_set.get_from_id(set_id, workgroup_name, workbook_name)
    serialised_set = services.reaction_set.to_dict([r_set])[0]

    return render_template(
        "reaction_set.html",
        reaction_set=serialised_set,
        number_of_reactions=len(r_set.reactions),
        workgroup_name=workgroup_name,
        workbook_name=workbook_name,
        reactor_dimensions=r_set.reactor_dimensions,
    )


@reaction_set_bp.route("/click_and_drag")
def click_and_drag():
    return render_template("click-and-drag.html")


@reaction_set_bp.route("/update_reaction_set", methods=["POST"])
def update_reaction_set() -> Response:
    """
    Returns serialised reaction set for storage on front end. Used after saving multiple wells to apply changes in
    front end
    """
    set_id = request.json.get("set_id")
    workgroup_name = request.json.get("workgroup")
    workbook_name = request.json.get("workbook")
    print(set_id, workgroup_name, workbook_name)

    r_set = services.reaction_set.get_from_id(set_id, workgroup_name, workbook_name)
    serialised_set = services.reaction_set.to_dict([r_set])[0]

    return serialised_set


@reaction_set_bp.route("/new_reaction_set", methods=["GET", "POST"])
@login_required
def new_reaction_set():
    data = request.get_json()
    workgroup_name = data.get("workgroup", None)
    workbook_name = data.get("workbook", None)
    reactor_dimensions = data.get("reactorDimensions", None)
    # probs need some checks here to make sure necessary data is present

    new_set_id = data.get("setID", None)
    new_set_name = data.get("setName")

    workgroup = services.workgroup.from_name(workgroup_name)
    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )
    creator = services.person.from_current_user_email()

    number_of_reactions = get_number_of_reactions(reactor_dimensions)

    # new_set_id = services.reaction_set.next_id_in_workbook(workbook.id)
    reaction_table = services.reaction.empty_reaction_table()
    summary_table = services.summary.empty_summary_table()

    reactions = []
    for i in range(number_of_reactions):
        name_and_id = new_set_id + "-" + str(i + 1)
        # reaction_id = services.reaction.get_next_reaction_id_for_workbook(workbook.id)
        reactions.append(
            services.reaction.add(
                # use same value for reaction_id and name
                name=name_and_id,
                reaction_id=name_and_id,
                creator=current_user,
                workbook_id=workbook.id,
                reaction_table=reaction_table,
                summary_table=summary_table,
                reaction_smiles="",
            )
        )

    services.reaction_set.add(
        name=new_set_name,
        set_id=new_set_id,
        creator=creator,
        workbook=workbook,
        workgroup=workgroup,
        reactions=reactions,
        reactor_dimensions=reactor_dimensions,
    )

    return jsonify({"feedback": "Set Created!"})


def get_number_of_reactions(reactor_dimensions):
    if reactor_dimensions.get("reactorType") == "carousel":
        return reactor_dimensions.get("numberOfReactions")
    elif reactor_dimensions.get("reactorType") == "well-plate":
        return reactor_dimensions.get("rows") * reactor_dimensions.get("columns")


@reaction_set_bp.route("/import_from_reactwise", methods=["POST", "GET"])
def import_from_reactwise():
    step_id = request.json.get("reactwiseID", None)
    workbook_name = request.json.get("workbook", None)
    workgroup_name = request.json.get("workgroup", None)

    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )

    data = ReactWiseStep(int(step_id))
    data.load_step_experiments()
    creator = services.user.person_from_current_user()

    set_name = "ReactWise Set " + str(step_id)

    # check if set already exists
    set_obj = services.reaction_set.get_from_names(
        set_name, workgroup_name, workbook_name
    )

    if not set_obj:
        reactions = []

        for reactwise_id, details in data.experimental_details.items():
            reaction_id = services.reaction.get_next_reaction_id_for_workbook(
                workbook.id
            )
            reaction = services.reaction.add(
                name="reactwise-" + reactwise_id,
                creator=creator,
                reaction_id=reaction_id,
                workbook_id=workbook.id,
                reaction_table={},
                summary_table={},
                reaction_smiles="OB(O)C1=CC=CC=C1.FC2=CC=C(Br)C=C2>>FC3=CC=C(C=C3)C4=CC=CC=C4",
            )
            reactions.append(reaction)

        set_id = services.reaction_set.next_id_in_workbook(workbook.id)
        services.reaction_set.add(
            name=set_name,
            set_id=set_id,
            creator=creator,
            workbook=workbook,
            reactions=reactions,
        )

    return jsonify(
        {
            "redirect_url": url_for(
                "reaction_set.reaction_set",
                set_name=set_name,
                workgroup_name=workgroup_name,
                workbook_name=workbook_name,
            )
        }
    )


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
