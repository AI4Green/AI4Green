import json
from datetime import datetime

import pytz
from compound_database import Compound_database_extraction as CDE
from rdkit import Chem
from sources import models, services
from sources.extensions import db


def add_test_user(
    username: str, email: str, fullname: str, password: str, verified: bool = False
):
    """
    Adds a test user to the database.
    Args:
        username (str): The username of the test user.
        email (str): The email address of the test user.
        fullname (str): The full name of the test user.
        password (str): The password of the test user.
        verified (bool): Whether the test user is verified.
    Returns:
        p (models.Person): The added user
    """
    # create and add person
    p = models.Person()
    db.session.add(p)

    # specify and add user
    services.user.add(username, email, fullname, password, p)

    if verified:
        user = services.user.from_email(email)
        user.is_verified = True

    db.session.commit()

    return p


def insert_test_data():
    time_of_creation = datetime.now(pytz.timezone("Europe/London")).replace(tzinfo=None)
    # make a person and a user
    CDE.seed_role_types()

    # add verified user
    # create and add person
    p1 = add_test_user(
        username="test_username",
        email="test_user@test.com",
        fullname="Gloria Testeban",
        password="test_pw",
        verified=True,
    )

    # add unverified user
    p2 = add_test_user(
        username="not_verified",
        email="not_verified@test.com",
        fullname="Test Daley",
        password="not_verified",
        verified=False,
    )

    p3 = add_test_user(
        username="password_reset",
        email="password_reset@test.com",
        fullname="TEST",
        password="password_reset",
        verified=True,
    )

    # Make an institution, workgroup and workbook
    institution1 = models.Institution.create(
        name="Test University", time_of_creation=time_of_creation
    )

    workgroup = models.WorkGroup.create(
        name="Test-Workgroup", institution=institution1.id, principal_investigator=[p1]
    )

    workbook = models.WorkBook.create(
        name="Test-Workbook", group=workgroup.id, users=[p1], abbreviation="TW1"
    )

    # add some compounds to the database
    compound1 = models.Compound(
        name="Testoic Acid",
        cas="123-45-5",
        cid=-10,
        hphrase="H900-H901",
        molec_weight=123,
        smiles="C",
    )
    compound2 = models.Compound(
        name="Testamine",
        cas="123-45-7",
        cid=-2,
        hphrase="H900",
        molec_weight=101,
        smiles="CC",
    )
    compound3 = models.Compound(
        name="Testide",
        cas="123-45-8",
        cid=-3,
        hphrase="H901",
        molec_weight=150,
        smiles="CCC",
    )
    compound4 = models.Compound(
        name="Testol",
        cas="123-45-9",
        cid=-4,
        hphrase="H901",
        molec_weight=175,
        smiles="CCCC",
    )
    compound5 = models.Compound(
        name="Testanol",
        cas="123-46-0",
        cid=-5,
        hphrase="H900",
        molec_weight=200,
        smiles="CCCCC",
    )
    compound6 = models.Compound(
        name="Testproductine",
        cas="123-46-1",
        cid=-6,
        hphrase="H900",
        molec_weight=225,
        smiles="CP",
    )
    compound7 = models.Compound(
        name="Testproductone",
        cas="123-46-2",
        cid=-7,
        hphrase="H900",
        molec_weight=255,
        smiles="CCP",
    )
    for compound in [
        compound1,
        compound2,
        compound3,
        compound4,
        compound5,
        compound6,
        compound7,
    ]:
        # add inchi using SMILES then add to db
        mol = Chem.MolFromSmiles(compound.smiles)
        compound.inchi = Chem.MolToInchi(mol)
        db.session.add(compound)

    # take compound5 and add it to solvent table
    models.Solvent.create(
        name="Testanol",
        flag=2,
        hazard="H900",
        compound=[compound5],
        time_of_creation=time_of_creation,
    )

    models.HazardCode.create(code="H900", phrase="Caution testing", category="")
    models.HazardCode.create(code="H901", phrase="May cause tests", category="")

    # data for saving a reaction
    reaction_table = json.dumps(reaction_table_json_contents())

    summary_table = json.dumps(summary_json_contents())

    models.Reaction.create(
        creator=p1.id,
        time_of_creation=time_of_creation,
        reaction_id="TW1-001",
        name="a test reaction",
        description="testing a reaction involving compounds from the test family",
        reaction_class="Amide bond formation",
        reactants=["C", "CC"],
        products=["CP", "CPP"],
        reagents=["CCC", "CCCCC"],
        solvent=[-5],
        reaction_smiles="C.CC>>CP.CPP",
        workbooks=workbook.id,
        reaction_table_data=reaction_table,
        summary_table_data=summary_table,
        complete="not complete",
        status="active",
        reaction_type="STANDARD",
    )

    # add a polymer reaction
    models.Reaction.create(
        creator=p1.id,
        time_of_creation=time_of_creation,
        reaction_id="TW1-002",
        name="a test polymer reaction",
        description="testing a reaction involving compounds from the test family",
        reaction_class="",
        reactants=["*C*", "CC"],
        products=["CP", "CPP"],
        reagents=["CCC", "CCCCC"],
        solvent=[-5],
        reaction_smiles="*C*.CC>>CP.CPP",
        workbooks=workbook.id,
        reaction_table_data=reaction_table,
        summary_table_data=summary_table,
        complete="not complete",
        status="active",
        reaction_type="POLYMER",
    )

    # add a novel compound
    inchi_1 = "InChI=1S/C21H15NO/c23-21(15-8-2-1-3-9-15)22-20-18-12-6-4-10-16(18)14-17-11-5-7-13-19(17)20/h1-14H,(H,22,23)"
    models.NovelCompound.create(
        name="Starting material",
        cas="",
        workbook=workbook.id,
        molec_formula="C21H15NO",
        hphrase="Unknown",
        density=None,
        concentration=None,
        molec_weight=297.12,
        smiles="O=C(NC1=C2C=CC=CC2=CC2=C1C=CC=C2)C1=CC=CC=C1",
        inchi=inchi_1,
        inchikey="MRGYEWYAVKQUEA-UHFFFAOYSA-N",
    )

    # add a polymer novel compound
    models.PolymerNovelCompound.create(
        name="Polymer Compound",
        workbook=workbook.id,
        molec_formula="C",
        hphrase="Unknown",
        density=None,
        concentration=None,
        molec_weight=12.01,
        smiles="*C*",
    )

    seed_limited_element_sustainability_data()

    print("Database populated with test data")


def seed_limited_element_sustainability_data():
    models.Element.create(
        name="Carbon",
        symbol="C",
        remaining_supply=500,
        colour="lime",
    ),
    models.Element.create(
        name="Nitrogen",
        symbol="N",
        remaining_supply=500,
        colour="lime",
    ),
    models.Element.create(
        name="Oxygen",
        symbol="O",
        remaining_supply=500,
        colour="lime",
    ),
    models.Element.create(
        name="Phosphorus",
        symbol="P",
        remaining_supply=500,
        colour="lime",
    )


def reaction_table_json_contents():
    return {
        "reaction_smiles": "C.CC>>CP.CCP",
        "reaction_name": "test reaction name",
        "reaction_description": "testing routes and services",
        "reaction_class": "Amide bond formation",
        # reactant data
        "reactant_ids": [-1, -2],
        "reactant_masses": ["123", "101"],
        "reactant_masses_raw": ["123", "101"],
        "reactant_amounts": ["1.00", "10.0"],
        "reactant_amounts_raw": ["1.00", "10.0"],
        "reactant_volumes": ["0.10", "0.02"],
        "reactant_densities": ["1.204", "0.995"],
        "reactant_volumes_raw": ["0.10", "0.02"],
        "reactant_equivalents": ["1", "10"],
        "reactant_concentrations": ["-", "-"],
        "reactant_physical_forms": ["2", "6"],
        "limiting_reactant_table_number": 1,
        "reagent_ids": [-3, -4],
        "reagent_names": ["Testide", "Testol"],
        "reagent_amounts": ["0.20", "2.60"],
        "reagent_amounts_raw": ["0.20", "2.60"],
        "reagent_equivalents": ["0.2", "2.6"],
        "reagent_molecular_weights": ["150", "175"],
        "reagent_densities": ["", "0.703"],
        "reagent_concentrations": ["", ""],
        "reagent_masses": ["3.21", "297"],
        "reagent_masses_raw": ["3.2054", "296.706"],
        "reagent_volumes": ["", "0.42"],
        "reagent_volumes_raw": ["", "0.422475"],
        "reagent_physical_forms": ["1", "1"],
        "reagent_hazards": ["H220", "H225-H302-H304-H312-H315-H319-H336-H400-H410"],
        "solvent_names": ["Testanol"],
        "solvent_concentrations": ["0.50", "1.00"],
        "solvent_ids": [-5],
        "solvent_volumes": ["2", "1"],
        "solvent_physical_forms": ["2", "6"],
        "solvent_hazards": ["H225-H302-H319-H371", "H225-H302-H319-H335-H351"],
        "product_ids": [-6, -7],
        "product_masses": ["93.0", "18.0"],
        "product_masses_raw": ["93.0", "18.0"],
        "product_amounts": ["1.00", "1.00"],
        "product_amounts_raw": ["1.00", "1.00"],
        "product_mass_units": "g",
        "product_physical_forms": ["7", "6"],
        "main_product_table_number": 7,
        "main_product": 1,
        "amount_units": "mmol",
        "mass_units": "mg",
        "volume_units": "mL",
        "solvent_volume_units": "mL",
        "product_amount_units": "mmol",
    }


def summary_json_contents():
    return {
        "real_product_mass": "90",
        "unreacted_reactant_mass": "2",
        "reaction_temperature": "70",
        "element_sustainability": "3",
        "batch_flow": "Batch",
        "isolation_method": "2",
        "catalyst_used": "Catalyst or enzyme",
        "catalyst_recovered": "Not recovered catalyst",
        "radio_buttons": ["hplc", "massSpec"],
        "custom_protocol1": "",
        "custom_protocol2": "",
        "other_hazards_text": "Weigh out reactants in fumehood",
        "researcher": "",
        "supervisor": "",
    }
