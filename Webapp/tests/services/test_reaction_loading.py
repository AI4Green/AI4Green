from flask import Flask
from flask.testing import FlaskClient
from sources.extensions import db
from test_reaction_table import make_url, assert_reaction_table_response_for_test_compounds
from pathlib import Path
from tests.utils import login
import pickle
from pytest_mock import MockFixture
from sources import services, models

def test_reload_reaction_v1_5(app):
    pass


def test_reload_reaction_v1_6(app: Flask, client: FlaskClient):
    login(client)
    with app.app_context():
        reaction_file = Path(__file__).resolve().parent.parent / "data" / "reaction_database_object_v1_6.pickle"
        with open(reaction_file, "rb") as f:
            reaction = pickle.load(f)
            db.session.add(reaction)
            db.session.commit()

            # test sketcher load
            sketcher_url = f"sketcher/{reaction.workbook.WorkGroup.name}/{reaction.workbook.name}/{reaction.reaction_id}/no"
            response = client.get(sketcher_url)
            assert response.status_code == 200

            # test reaction table: use only one product as reaction.products[-1] is not in database
            table_url = make_url(reactants=",".join(reaction.reactants), products=reaction.products[0] + "," + reaction.products[0])
            response = client.get(table_url)
            assert_reaction_table_response_for_test_compounds(response)


def test_clone_reaction(app: Flask, client: FlaskClient, mocker: MockFixture):
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
        assert services.reaction.get_from_reaction_id_and_workbook_id("TW1-002", 1) is not None


def test_delete_reaction(client: FlaskClient, app: Flask):
    login(client)
    with app.app_context():
        reaction = services.reaction.get_from_reaction_id_and_workbook_id("TW1-001", 1)
        response = client.get(f"delete_reaction/{reaction.reaction_id}/{reaction.workbook.WorkGroup.name}/{reaction.workbook.name}")
        assert response.status_code == 302
        assert reaction.status == "inactive"


def save_reaction_version(app: Flask):
    """
    function to save a test Reaction from the database. Reactions from each major version of
    AI4Green should be saved and tested for reload to ensure old reactions are compatible with
    new updates.
    """
    with app.app_context():
        version = "CHANGE_ME"
        reaction = (
        db.session.query(models.Reaction)
        .filter(models.Reaction.name == "a test reaction")
        .join(models.WorkBook)
        .filter(models.WorkBook.id == 1)
        .first()
        )

        p = pickle.dumps(reaction)
        with open("data/reaction_database_object" + version + ".pickle", "wb+") as f:
            f.write(p)

