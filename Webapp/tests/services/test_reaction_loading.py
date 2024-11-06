from flask import Flask
from flask.testing import FlaskClient
from sources import services, models
from sources.extensions import db
from pathlib import Path
from tests.utils import login
import pickle

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

            sketcher_url = f"sketcher/{reaction.workbook.WorkGroup.name}/{reaction.workbook.name}/{reaction.reaction_id}/no"
            response = client.get(sketcher_url)
            assert response.status_code == 200


def test_clone_reaction():
    pass

def test_delete_reaction():
    pass

def test_save_reaction_version(client, app):
    """ function to save Reaction from a specified version of the database."""
    with app.app_context():
        reaction = (
        db.session.query(models.Reaction)
        .filter(models.Reaction.name == "a test reaction")
        .join(models.WorkBook)
        .filter(models.WorkBook.id == 1)
        .first()
        )

        p = pickle.dumps(reaction)

        with open("data/reaction_database_object_v1_6.pickle", "wb+") as f:
            f.write(p)

