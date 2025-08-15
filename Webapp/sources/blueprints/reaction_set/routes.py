from typing import Dict, List, Union

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
    step_id = request.json.get("stepID", None)
    step_name = request.json.get("stepName", None)
    workbook_name = request.json.get("workbook", None)
    workgroup_name = request.json.get("workgroup", None)

    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )

    rw_step = services.reactwise.ReactWiseStep(int(step_id))
    rw_step.load_step_experiments()
    creator = services.person.from_current_user_email()

    set_name = step_name

    # check if set already exists
    set_obj = services.reaction_set.get_from_names(
        set_name, workgroup_name, workbook_name
    )

    if not set_obj:
        reactions = []

        for reactwise_id, details in rw_step.experimental_details.items():
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
