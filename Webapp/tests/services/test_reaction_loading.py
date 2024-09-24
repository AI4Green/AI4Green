from sources import services, models
from sources.extensions import db
import pickle

def test_reload_reaction_v1_5():
    pass

def test_reload_reaction_v1_6():
    pass

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

