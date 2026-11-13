"""
Extra seed data for the end-to-end tests, added on top of tests/populate_database.py.

Production databases get the COSHH input types and the H&S Admin role from an
Alembic migration, but the test database is built with ``db.create_all()``, so
they are added here. Every function must run inside an app context.
"""

from datetime import datetime

from e2e.seed_data import (
    COSHH_SECTION,
    COSHH_TEMPLATE,
    COSHH_TEXT_FIELD,
    HS_ADMIN_USERNAME,
    INVENTORY_LOCATION,
    INVENTORY_OTHER_LOCATION,
    INVENTORY_WORKGROUP,
    ROLE_USER_PASSWORD,
    STOCKED_COMPOUND_CAS,
    TEST_USERNAME,
)
from sources import models, services
from sources.extensions import db
from tests import populate_database

HS_ADMIN_ROLE = "H&S Admin"

# copied from migrations/versions/13d40f5a5809_coshh_template.py, ids included
COSHH_INPUT_TYPES = [
    "content",
    "text",
    "description",
    "number",
    "reaction",
    "multi-reaction",
    "radio",
    "multiple",
    "file",
    "image-file",
    "header",
    "chemicaldisposal",
    "sortable-list",
    "projectgroup-plan",
    "projectgroup-hazard",
    "yield",
    "multi-yield",
    "greenmetrics",
    "multi-greenmetrics",
    "datetime",
    "formatted-text-input",
]


def _user(username: str) -> models.User:
    """
    Look up a seeded user by username.

    Args:
        username: The user's username.

    Returns:
        The matching user.
    """
    return db.session.query(models.User).filter_by(username=username).one()


def _workgroup(name: str) -> models.WorkGroup:
    """
    Look up a seeded workgroup by name.

    Args:
        name: The workgroup's name.

    Returns:
        The matching workgroup.
    """
    return db.session.query(models.WorkGroup).filter_by(name=name).one()


def add_coshh_data() -> None:
    """
    Add what the COSHH tests need: the input types, an H&S admin who can design
    templates, and a template with one section holding one optional text field.
    """
    for input_type_id, title in enumerate(COSHH_INPUT_TYPES, start=1):
        db.session.add(models.InputType(id=input_type_id, title=title))
    hs_admin_role = models.Role.create(
        name=HS_ADMIN_ROLE,
        role_description="H&S Admin with permissions to design COSHH templates",
    )

    populate_database.add_test_user(
        username=HS_ADMIN_USERNAME,
        email=f"{HS_ADMIN_USERNAME}@test.com",
        fullname="E2E Health and Safety Admin",
        password=ROLE_USER_PASSWORD,
        verified=True,
    )
    hs_admin = _user(HS_ADMIN_USERNAME)
    hs_admin.role = hs_admin_role.id

    template = models.COSHHTemplate(
        name=COSHH_TEMPLATE,
        description="Template seeded for the end-to-end tests",
        time_of_creation=datetime.now(),
        creator_id=hs_admin.id,
        institution_id=db.session.query(models.Institution).first().id,
    )
    section = models.Section(name=COSHH_SECTION, sort_order=1, template=template)
    db.session.add(
        models.Field(
            name=COSHH_TEXT_FIELD,
            sort_order=1,
            mandatory=False,
            section=section,
            input_type_id=COSHH_INPUT_TYPES.index("text") + 1,
        )
    )
    db.session.commit()


def add_inventory_data() -> None:
    """
    Give the inventory workgroup two storage locations and stock it with one
    compound, so tests can act on containers without adding them first.
    """
    workgroup = _workgroup(INVENTORY_WORKGROUP)
    for name in (INVENTORY_LOCATION, INVENTORY_OTHER_LOCATION):
        db.session.add(models.Location(workgroup_id=workgroup.id, name=name))

    compound = (
        db.session.query(models.Compound).filter_by(cas=STOCKED_COMPOUND_CAS).one()
    )
    db.session.add(
        models.WorkgroupCompound(
            workgroup_id=workgroup.id, compound_id=compound.id, threshold_unit="mL"
        )
    )
    db.session.commit()


def add_inventory_container(lot_number: str, quantity: int = 500) -> None:
    """
    Add a container of the stocked compound to the first inventory location,
    recorded in the ledger as received by the main test user.

    Args:
        lot_number: Lot number for the container, which tests use to find it.
        quantity: Amount in the container, in mL.
    """
    workgroup = _workgroup(INVENTORY_WORKGROUP)
    location = (
        db.session.query(models.Location)
        .filter_by(workgroup_id=workgroup.id, name=INVENTORY_LOCATION)
        .one()
    )
    workgroup_compound = (
        db.session.query(models.WorkgroupCompound)
        .join(models.Compound)
        .filter(
            models.WorkgroupCompound.workgroup_id == workgroup.id,
            models.Compound.cas == STOCKED_COMPOUND_CAS,
        )
        .one()
    )
    services.inventory_ledger.create_inventory_item(
        workgroup_compound.id,
        location.id,
        quantity,
        "mL",
        _user(TEST_USERNAME).id,
        lot_number=lot_number,
        supplier="E2E Supplies",
    )
    db.session.commit()
