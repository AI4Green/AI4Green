import json
from typing import Dict, List, Union

from flask import (
    Response,
    current_app,
    jsonify,
    redirect,
    render_template,
    request,
    session,
    url_for,
)
from flask_login import current_user, login_required
from sources import models, services
from sources.extensions import db

from . import reaction_set_bp


@reaction_set_bp.route("/reaction_set/<workgroup_name>/<workbook_name>/<set_name>")
@login_required
def reaction_set(workgroup_name, workbook_name, set_name):
    r_set = services.reaction_set.from_names(set_name, workgroup_name, workbook_name)
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


def rw_reactor_dimensions(number_of_reactions):
    if number_of_reactions <= 12:  # less than 12 reactions
        return {"reactorType": "carousel", "numberOfReactions": number_of_reactions}
    else:
        wellplates = {
            24: {"rows": 4, "cols": 6},
            48: {"rows": 6, "cols": 8},
            96: {"rows": 8, "cols": 12},
        }
        # Pick the smallest standard plate that fits
        for size, dims in wellplates.items():
            if number_of_reactions <= size:
                return {
                    "reactorType": "well plate",
                    "rows": dims["rows"],
                    "columns": dims["cols"],
                }


@reaction_set_bp.route("/import_from_reactwise", methods=["POST", "GET"])
def import_from_reactwise():
    step_id = request.json.get("stepID", None)
    step_name = request.json.get("stepName", None)
    workbook_name = request.json.get("workbook", None)
    workgroup_name = request.json.get("workgroup", None)
    workgroup = services.workgroup.from_name(workgroup_name)
    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )

    rw_step = services.reactwise.ReactWiseStep(int(step_id))
    creator = services.person.from_current_user_email()

    # check if set already exists
    set_obj = services.reaction_set.from_names(step_name, workgroup_name, workbook_name)
    # handle error here lmao

    if not set_obj:
        set_id = services.reaction_set.next_id_in_workbook(workbook.id)
        reactions = []
        unknown_solvents = {}
        unknown_fields = []
        novel_compounds = []
        reaction_components = []

        # these fields are added to the reaction table dict
        known_fields = ["solvent", "time", "temperature"]  # any more to add?

        for reactwise_experiment_id, details in rw_step.experimental_details.items():
            # fix me
            unknown_keys = [
                (param, y.get("unit", "") if isinstance(y, dict) else "")
                for param, y in details.items()
                if param.lower() not in known_fields
            ]
            unknown_fields.extend(unknown_keys)

            # TODO: move exp details logic all to assign variables route?

            # load new reaction table to load into reaction
            new_reaction_table_data = services.reaction_table.new()
            reaction_id = services.reaction.get_next_reaction_id_for_workbook(
                workbook.id
            )

            reaction = services.reaction.add(
                name="reactwise-" + reactwise_experiment_id,
                creator=creator,
                reaction_id=reaction_id,
                workbook_id=workbook.id,
                reaction_table=new_reaction_table_data,
                summary_table={},
                reaction_smiles=rw_step.reaction_smiles,
            )

            # gets solvents
            solvent = services.reactwise.extract_reaction_solvents(
                reaction.reaction_id, details, unknown_solvents
            )

            # add reactants/products/solvents
            reaction.reactants = rw_step.reactants
            reaction.products = rw_step.products
            reaction.solvent = [solvent.id if solvent else ""]

            reaction_components.extend(rw_step.reactants)
            reaction_components.extend(rw_step.products)

            # now process the reaction data using the reaction table class
            reaction_table = services.reaction_table.ReactionTable(
                reaction, workgroup, workbook, "no", "no"
            )
            # this step updates the products/reagents from reaction smiles in the reaction table dict
            # if there is a novel compound, they are added as blank strings and the novel_compound responses are stored
            reaction_table.update(rw_step.reaction_smiles, {})
            reaction_table.update_reaction_table_data()
            novel_compounds.extend(reaction_table.novel_compounds)

            # reaction_table_data = reaction_table.reaction_table_data

            # this step loads reaction temp/time to the reaction table dict
            services.reactwise.extract_temp_time_fields(
                details, reaction_table.reaction_table_data
            )
            print(type(reaction.reaction_table_data), "1")

            reaction_table.save()
            print(type(reaction.reaction_table_data), "2")

            reactions.append(reaction)

        reactor_dimensions = rw_reactor_dimensions(len(reactions))

        services.reaction_set.add(
            name=step_name,
            set_id=set_id,
            creator=creator,
            workgroup=workgroup,
            workbook=workbook,
            reactions=reactions,
            reactor_dimensions=reactor_dimensions,
        )
        if unknown_solvents or unknown_fields:
            # store as session variables for now, there must be a better way
            session["reactwise_import"] = {
                "unknown_fields": list(set(unknown_fields)),
                "unknown_solvents": unknown_solvents,
                "set_id": set_id,
                "workgroup_name": workgroup_name,
                "workbook_name": workbook_name,
                "novel_compounds": list(set(novel_compounds)),
                "reaction_smiles": rw_step.reaction_smiles,
                "reaction_components": list(set(reaction_components)),
                "experimental_details": rw_step.experimental_details,
                "set_name": step_name,
            }
            print("SETUP", rw_step.experimental_details)

            return jsonify(
                {
                    "feedback": "unknown fields",
                }
            )
        return jsonify({"feedback": "Set Created!"})
    else:
        return jsonify({"feedback": "Set Created!"})


@reaction_set_bp.route("/assign_reactwise_fields", methods=["GET", "POST"])
def assign_reactwise_fields():
    # requires session variable, keeps in session in case of reload
    rw_import = session.get("reactwise_import", {})
    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        rw_import["workgroup_name"], rw_import["workbook_name"]
    )
    novel_compounds = []

    for smiles in rw_import["novel_compounds"]:
        novel_compounds.append(
            {
                "name": services.compound.iupac_convert(smiles),
                "mw": services.all_compounds.mol_weight_from_smiles(smiles),
                "smiles": smiles,
                "image_base64": services.utils.smiles_to_base64(smiles),
            }
        )

    reaction_image = services.utils.reaction_to_base64(rw_import["reaction_smiles"])

    sol_rows = services.solvent.get_workbook_list(workbook)

    return render_template(
        "reactions/assign_reactwise_fields.html",
        unknown_variables=rw_import["unknown_fields"],
        unknown_solvents=rw_import["unknown_solvents"],
        number_of_solvents=len(rw_import["unknown_solvents"]),
        number_of_variables=len(rw_import["unknown_fields"]),
        sol_rows=sol_rows,
        workgroup_name=rw_import["workgroup_name"],  # fix hidden input bug
        workbook_name=rw_import["workbook_name"],
        novel_compounds=novel_compounds,
        reaction_scheme_image=reaction_image,
        reaction_components=rw_import["reaction_components"],
        reaction_set_id=rw_import["set_id"],
        set_name=rw_import["set_name"],
        experimental_details=rw_import["experimental_details"],
        reaciton_smiles=rw_import["reaction_smiles"],
    )


# possibly patch
@reaction_set_bp.route("/assign_reactwise_variables", methods=["POST"])
def assign_reactwise_variables():
    assignments = request.json.get("assignments")
    workbook_name = request.json.get("workbook_name")
    workgroup_name = request.json.get("workgroup_name")
    reaction_set_id = request.json.get("reaction_set_id")
    exp_details = request.json.get("exp_details")

    rset = services.reaction_set.get_from_id(
        reaction_set_id, workgroup_name, workbook_name
    )

    success = []
    fail = []

    for reaction_number, details in exp_details.items():
        for variable in details:
            if variable in assignments:
                ls = assignments.get(variable)
                print("ls", ls)
                reaction = [
                    r
                    for r in rset.reactions
                    if r.name == "reactwise-" + reaction_number
                ][0]
                reaction_table_data = json.loads(reaction.reaction_table_data)

                component_idx = reaction.reactants.index(ls["component"])
                value = float(details[variable].get("value"))

                if ls.get("type") == "mole-fraction":
                    scaled_value = value / 100
                    reaction_table_data["reactant_equivalents"][
                        component_idx
                    ] = scaled_value
                    success.append(variable)
                    reaction.reaction_table_data = json.dumps(reaction_table_data)

                print(type(reaction.reaction_table_data), "4")

                # include extra logic for other types

    db.session.commit()

    return jsonify({"success": success, "failed": fail})
