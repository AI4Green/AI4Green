from flask import Response, jsonify, request
from flask_login import login_required
from sources import services
from sources.auxiliary import abort_if_user_not_in_workbook

from . import novel_compound_bp


@novel_compound_bp.route("/_novel_compound", methods=["GET", "POST"])
@login_required
def novel_compound() -> Response:
    """
    Adds a novel compound to the database and links it to the workbook of the current reaction.
    The novel compound data is saved from the data the user enters in the form.
    This function is triggered either from the 'add new reagent/solvent to database'
    or when a novel structure is drawn in the sketcher.

    Returns:
        Flask Response with feedback about the operation.
    """
    # get the active workbook and verify user belongs
    workgroup_name = str(request.form["workgroup"])
    workbook_name = str(request.form["workbook"])
    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )
    abort_if_user_not_in_workbook(workgroup_name, workbook_name, workbook)

    # check all the data the user has supplied to ensure it is valid and no unique constraints are broken.
    new_compound = services.novel_compound.NewNovelCompound(workbook)
    validation_results = new_compound.validate()
    if validation_results == "validation_successful":
        # create novel compound
        nc = services.novel_compound.add(
            new_compound.name,
            new_compound.cas,
            new_compound.mol_formula,
            new_compound.mol_weight,
            new_compound.density,
            new_compound.concentration,
            new_compound.hazard_codes,
            new_compound.smiles,
            new_compound.inchi,
            new_compound.inchi_key,
            new_compound.workbook.id,
        )

        # if the user has added the compound as a novel solvent, we add this to the solvent table too
        component_type = request.form["component"]
        if component_type == "solvent":
            services.solvent.add(new_compound.name, new_compound.hazard_codes, nc)

    return jsonify({"feedback": new_compound.feedback})
