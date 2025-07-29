"""
This module receives a reaction from Marvin JS as a
GET request and renders the reaction table template
"""

from flask import json, jsonify, render_template, request
from flask_login import login_required
from sources import models, services
from sources.dto import ReactionNoteSchema

from . import reaction_table_bp


@reaction_table_bp.route("/autoupdate_reaction_table", methods=["GET", "POST"])
# @workbook_member_required
def autoupdate_reaction_table():
    """
    I guess the idea is to take any smiles from the front end and generate the reaction table

    Needs to get reaction from db too, but ignore for now
    """
    # get user workbook
    demo = request.json.get("demo")
    tutorial = request.json.get("tutorial")
    workbook = None
    reaction = None
    polymer_indices = []
    if demo != "demo" and tutorial != "yes":
        workgroup = request.json.get("workgroup")
        workbook_name = request.json.get("workbook")
        workbook = services.workbook.get_workbook_from_group_book_name_combination(
            workgroup, workbook_name
        )
        reaction_id = request.json.get("reaction_id")
        reaction = services.reaction.get_from_reaction_id_and_workbook_id(
            reaction_id, workbook.id
        )

    reaction_smiles = request.json.get("reaction_smiles")
    if not reaction_smiles:
        return jsonify({"error": "Missing data!"})

    polymer_indices = request.json.get("polymer_indices")

    reaction_table = services.reaction_Table.ReactionTable(reaction.reaction_table_data)

    print(reaction.reaction_table_data)

    (
        reactants_smiles_list,
        product_smiles_list,
    ) = services.reaction_table.get_reactants_and_products_list(reaction_smiles)
    reactants = [
        services.compound.SketcherCompound(
            smiles=x,
            idx=idx + 1,
            polymer_indices=polymer_indices,
            workbook=workbook,
            demo=demo,
            reaction_component="Reactant",
            reaction_component_idx=idx + 1,
        )
        for idx, x in enumerate(reactants_smiles_list)
    ]
    number_of_reactants = len(reactants)

    # add for reagent support
    # reagents = [
    #     SketcherCompound(
    #         smiles=x,
    #         idx=idx + 1 + number_of_reactants,
    #         polymer_indices=polymer_indices,
    #         workbook=workbook,
    #         demo=demo,
    #         reaction_component="Reagent"
    #     )
    #     for idx, x in enumerate(reactants_smiles_list)
    # ]
    # number_of_reagents = len(reactants)

    products = [
        services.compound.SketcherCompound(
            smiles=y,
            idx=idx + 1 + number_of_reactants,
            polymer_indices=polymer_indices,
            workbook=workbook,
            demo=demo,
            reaction_component="Product",
            reaction_component_idx=idx + 1,
        )
        for idx, y in enumerate(product_smiles_list)
    ]
    number_of_products = len(products)

    # check errors first add reagents here for reagent support
    for compound_group in (reactants, products):
        # now handles co polymer, dummy atom and invalid molecule errors
        error = services.compound.check_compound_errors(compound_group)
        if error:
            return error

        novel_compound_row = services.compound.check_novel_compounds(compound_group)
        if novel_compound_row:
            # this should be changed i think
            return novel_compound_row

    identifiers = []

    # Solvents - keep solvents that are not novel compounds or are novel compounds within the current workbook
    if demo == "demo":
        sol_rows = services.solvent.get_default_list()
    else:
        sol_rows = services.solvent.get_workbook_list(workbook)

    r_class = None

    if not polymer_indices:
        polymer_indices = list()
        r_class = services.reaction_classification.classify_reaction(
            reactants_smiles_list, product_smiles_list
        )

    reaction_table_html = "reactions/_reaction_table.html"

    # Now it renders the reaction table template
    reaction_table = render_template(
        reaction_table_html,
        reactants=reactants,
        # reagents=reagents,
        number_of_reactants=number_of_reactants,
        number_of_products=number_of_products,
        number_of_reagents=0,
        identifiers=identifiers,
        reactant_table_numbers=[],
        products=products,
        # units=default_units,
        # product_intended_dps=product_data["intended_dps"],
        reagent_table_numbers=[],
        reaction_table_data="",
        summary_table_data="",
        sol_rows=sol_rows,
        reaction=reaction,
        reaction_class=r_class,
        polymer_indices=polymer_indices,
    )
    return jsonify({"reactionTable": reaction_table})


@reaction_table_bp.route("/reload_reaction_table", methods=["GET", "POST"])
def reload_reaction_table():
    # load variables from front end
    workbook = request.json.get("workbook")
    workgroup = request.json.get("workgroup")
    reaction_id = request.json.get("reaction_id")
    demo = request.json.get("demo")
    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup, workbook
    )
    reaction = services.reaction.get_from_reaction_id_and_workbook_id(
        reaction_id, workbook.id
    )

    reaction_table = services.reaction_table.ReactionTable(
        reaction, workgroup, workbook, demo, "no"
    )
    print(reaction_table.reactants, reaction_table.products)

    return reaction_table.render()


@reaction_table_bp.route("/_save_reaction_note", methods=["POST"])
@login_required
@reaction_table_bp.doc(security="sessionAuth")
def save_reaction_note():
    """
    Saves an reaction_note to the reaction object

    Returns:
        flask.Response: A JSON response with the reaction_note object
    """
    workgroup_name = request.form["workgroup"]
    workbook_name = request.form["workbook"]
    workbook_object = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )
    reaction_id = request.form["reactionID"]
    reaction = services.reaction.get_from_reaction_id_and_workbook_id(
        reaction_id, workbook_object.id
    )
    reaction_note_text = request.form["reactionNoteText"]
    author = services.person.from_current_user_email()
    new_addendum = services.reaction.add_addendum(reaction, reaction_note_text, author)
    schema = ReactionNoteSchema()
    return jsonify({"reaction_note": schema.dump(new_addendum)})
