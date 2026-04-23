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
    input_type_table = sa.table(
        "InputType",
        sa.column("id", sa.Integer),
        sa.column("title", sa.String),
    )

    op.bulk_insert(
        input_type_table,
        [
            {"id": 1, "title": "Content"},
            {"id": 2, "title": "Text"},
            {"id": 3, "title": "Description"},
            {"id": 4, "title": "Number"},
            {"id": 5, "title": "Reaction Scheme"},
            {"id": 6, "title": "Multi Reaction Scheme"},
            {"id": 7, "title": "Radio"},
            {"id": 8, "title": "Multiple"},
            {"id": 9, "title": "File"},
            {"id": 10, "title": "Image File"},
            {"id": 11, "title": "Header"},
            {"id": 12, "title": "Chemical Disposal Table"},
            {"id": 13, "title": "Sortable List"},
            {"id": 14, "title": "Project Group Plan Table"},
            {"id": 15, "title": "Project Group Hazard Table"},
            {"id": 16, "title": "Yield Table"},
            {"id": 17, "title": "Multi Yield Table"},
            {"id": 18, "title": "Green Metrics Table"},
            {"id": 19, "title": "Multi Green Metrics Table"},
            {"id": 20, "title": "Date and Time"},
            {"id": 21, "title": "Formatted Text Input"},
        ],
    )


def downgrade_():
    # Clean up the seeded data if we roll back
    op.execute('DELETE FROM "InputType" WHERE id BETWEEN 1 AND 21')


def upgrade_audit_log():
    pass


def downgrade_audit_log():
    pass
