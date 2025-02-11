from datetime import datetime
from typing import Tuple

from flask import Flask
from flask.testing import FlaskClient
from flask_login import current_user
from sources import models, services
from sources.extensions import db
from tests.utils import login


def test_check_reaction_for_controlled_substances(app: Flask, client: FlaskClient):
    """
    Tests whether controlled substances are correctly identified in a given reaction
    """

    reaction_name, workbook_id = create_controlled_substance_reaction(app)
    with app.app_context():
        with client:
            login(client)
            reaction = services.reaction.get_from_name_and_workbook_id(
                reaction_name, workbook_id
            )
            checks = (
                services.controlled_substances.check_reaction_for_controlled_substances(
                    reaction
                )
            )

            # check return value
            return_value = [
                ["InChI=1S/C4H8Cl2S/c5-1-3-7-4-2-6/h1-4H2"],
                [],
                [
                    "InChI=1S/C6H6N12O12/c19-13(20)7-1-2-8(14(21)22)5(7)6-9(15(23)24)3(11(1)17(27)28)4(10(6)16(25)26)12(2)18(29)30/h1-6H"
                ],
                [],
            ]
            assert checks == return_value

            # check usages were recorded
            entry = services.controlled_substances.list_all()
            assert len(entry) == 2

            # check duplicates are correctly handled
            duplicate_checks = (
                services.controlled_substances.check_reaction_for_controlled_substances(
                    reaction
                )
            )
            assert duplicate_checks is None

            db.session.delete(reaction)
            db.session.commit()


def create_controlled_substance_reaction(app: Flask) -> Tuple[str, str]:
    """
    Creates a reaction containing controlled substances and adds to the database.
    Returns:
        reaction name and workbook id, for finding the reaction in the database
    """
    with app.app_context():
        workbook = services.workbook.get_workbook_from_group_book_name_combination(
            "Test-Workgroup", "Test-Workbook"
        )
        reaction = models.Reaction.create(
            creator=1,
            time_of_creation=datetime.now(),
            reaction_id="TW1-404",
            name="Chemical Weapon Test",
            description="testing a reaction involving chemical weapons",
            reaction_class="Amide bond formation",
            reactants=["ClCCSCCCl"],  # in db
            products=[
                "[O-][N+](N1C2N([N+]([O-])=O)C3N([N+]([O-])=O)C(N(C2N4[N+]([O-])=O)[N+]([O-])=O)C4N([N+]([O-])=O)C13)=O"
            ],  # novel compound
            reagents=[""],
            solvent=[-5],
            reaction_smiles="ClCCSCCCl>>O-][N+](N1C2N([N+]([O-])=O)C3N([N+]([O-])=O)C(N(C2N4[N+]([O-])=O)[N+]([O-])=O)C4N([N+]([O-])=O)C13)=O",
            workbooks=workbook.id,
            reaction_table_data="{}",
            summary_table_data="{}",
            complete="not complete",
            status="active",
        )
        db.session.add(reaction)
        db.session.commit()

        return reaction.name, workbook.id
