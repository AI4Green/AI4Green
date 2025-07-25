from flask import Response, current_app, jsonify, render_template
from flask_login import current_user, login_required
from sources import models, services
from sources.auxiliary import get_notification_number, get_workgroups
from sources.decorators import workbook_member_required
from sources.extensions import db

from . import reaction_constructor_bp  # imports the blueprint of the main route


@reaction_constructor_bp.route("/get_marvinjs_key", methods=["POST"])
def get_marvinjs_key():
    return jsonify({"marvinjs_key": current_app.config["MARVIN_JS_API_KEY"]})


# Go to the sketcher
@reaction_constructor_bp.route(
    "/sketcher/<workgroup>/<workbook>/<reaction_id>/<tutorial>", methods=["GET", "POST"]
)
@login_required
@workbook_member_required
def sketcher(
    workgroup: str, workbook: str, reaction_id: str, tutorial: str
) -> Response:
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    workbook_object = (
        db.session.query(models.WorkBook)
        .filter(models.WorkBook.name == workbook)
        .filter(models.WorkGroup.name == workgroup)
        .first()
    )
    reaction = (
        db.session.query(models.Reaction)
        .filter(models.Reaction.reaction_id == reaction_id)
        .filter(models.WorkBook.id == workbook_object.id)
        .first()
    )
    addenda = (
        db.session.query(models.ReactionNote)
        .join(models.Reaction)
        .filter(models.Reaction.id == reaction.id)
        .all()
    )

    if reaction.reaction_smiles:
        load_status = "loading"
    else:
        load_status = "loaded"

    return render_template(
        "reactions/reaction_constructor.html",
        reaction=reaction,
        load_status=load_status,
        demo="not demo",
        workgroups=workgroups,
        notification_number=notification_number,
        active_workgroup=workgroup,
        active_workbook=workbook,
        tutorial=tutorial,
        addenda=addenda,
    )


# Go to the sketcher tutorial
@reaction_constructor_bp.route("/sketcher_tutorial/<tutorial>", methods=["GET", "POST"])
def sketcher_tutorial(tutorial: str) -> Response:
    workgroups = []
    notification_number = 0
    if current_user.is_authenticated:
        workgroups = get_workgroups()
        notification_number = get_notification_number()
    return render_template(
        "reactions/reaction_constructor.html",
        reaction={
            "name": "Tutorial Reaction",
            "reaction_id": "TUT-001",
            "reaction_type": "STANDARD",
        },
        demo="not demo",
        workgroups=workgroups,
        notification_number=notification_number,
        active_workgroup=None,
        active_workbook=None,
        tutorial=tutorial,
    )
