"""populate_new_tables

Revision ID: 4f171d77ec98
Revises: 9445024a3863
Create Date: 2026-04-21 14:29:07.991914

"""
import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision = "4f171d77ec98"
down_revision = "9445024a3863"
branch_labels = None
depends_on = None


def upgrade(engine_name):
    globals()["upgrade_%s" % engine_name]()


def downgrade(engine_name):
    globals()["downgrade_%s" % engine_name]()


def upgrade_():
    # Define table structure locally for the migration
    section_type_table = sa.table(
        "SectionType", sa.column("id", sa.Integer), sa.column("name", sa.String)
    )

    # Data to seed
    op.bulk_insert(
        section_type_table,
        [
            {"id": 1, "name": "Activity Information"},
            {"id": 2, "name": "Hazard Identifications"},
            {"id": 3, "name": "SDS Reference"},
            {"id": 4, "name": "Safe Storage"},
            {"id": 5, "name": "Workplace Exposure Limits"},
            {"id": 6, "name": "Control Measures"},
            {"id": 7, "name": "Exposure Monitoring and Health Surveillance"},
            {"id": 8, "name": "Personal Protective Equipment (PPE)"},
            {"id": 9, "name": "Emergency Requirements"},
            {"id": 10, "name": "Additional Control Measures"},
        ],
    )


def downgrade_():
    # Clean up the seeded data if we roll back
    op.execute('DELETE FROM "SectionType" WHERE id BETWEEN 1 AND 10')


def upgrade_audit_log():
    pass


def downgrade_audit_log():
    pass
