from flask import Flask
from flask.testing import FlaskClient
from sources.extensions import db
from test_reaction_table import make_url, assert_reaction_table_response_for_test_compounds
from pathlib import Path
from tests.utils import login
import pickle
from pytest_mock import MockFixture
from sources import services, models


def test_clone_reaction(app: Flask, client: FlaskClient, mocker: MockFixture):
    """Tests clone_reaction functionality. Cloned reaction is removed from database at end of test."""
    login(client)
    with app.app_context():
        reaction = services.reaction.get_from_reaction_id_and_workbook_id("TW1-001", 1)
        mock_get_current_from_request_form = mocker.patch("sources.services.reaction.get_current_from_request_form")
        mock_get_current_from_request_form.return_value = reaction

        response = client.get("/clone_reaction", data={
            "reactionName": "Test Clone",
            "workbook": reaction.workbook.name,
            "workgroup": reaction.workbook.WorkGroup.name,
            "newReactionID": "TW1-002"
        })
        assert response.status_code == 200
        assert response.json["feedback"] == "New reaction made"
        new_reaction = services.reaction.get_from_reaction_id_and_workbook_id("TW1-002", 1)
        assert new_reaction is not None

        # delete cloned reaction for future tests
        db.session.delete(new_reaction)
        db.session.commit()


def test_delete_reaction(client: FlaskClient, app: Flask):
    """Tests delete_reaction functionality. Deleted reaction is labelled as active again at end of test."""
    login(client)
    with app.app_context():
        reaction = services.reaction.get_from_reaction_id_and_workbook_id("TW1-001", 1)
        response = client.get(f"delete_reaction/{reaction.reaction_id}/{reaction.workbook.WorkGroup.name}/{reaction.workbook.name}")
        assert response.status_code == 302
        assert reaction.status == "inactive"

        # undo delete for later tests
        reaction.status = "active"
        db.session.commit()


def save_reaction_version(app: Flask):
    """
    Function to save a test Reaction from the database. Reactions from each major version of
    AI4Green should be saved and tested for reload to ensure old reactions are compatible with
    new updates. Included here so that future versions can be saved in a consistent way.
    """
    with app.app_context():
        version = "vCHANGE_ME"
        reaction = (
        db.session.query(models.Reaction)
        .filter(models.Reaction.name == "Reloaded Reaction")
        .join(models.WorkBook)
        .filter(models.WorkBook.id == 1)
        .first()
        )
        db.session.expunge(reaction)

        p = pickle.dumps(reaction)
        with open("data/reaction_database_object_" + version + ".pickle", "wb+") as f:
            f.write(p)

