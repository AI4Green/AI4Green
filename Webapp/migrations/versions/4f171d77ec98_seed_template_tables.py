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
            {"id": 1, "title": "content"},
            {"id": 2, "title": "text"},
            {"id": 3, "title": "description"},
            {"id": 4, "title": "number"},
            {"id": 5, "title": "reaction"},
            {"id": 6, "title": "multi-reaction"},
            {"id": 7, "title": "radio"},
            {"id": 8, "title": "multiple"},
            {"id": 9, "title": "file"},
            {"id": 10, "title": "image-file"},
            {"id": 11, "title": "header"},
            {"id": 12, "title": "chemicaldisposal"},
            {"id": 13, "title": "sortable-list"},
            {"id": 14, "title": "projectgroup-plan"},
            {"id": 15, "title": "projectgroup-hazard"},
            {"id": 16, "title": "yield"},
            {"id": 17, "title": "multi-yield"},
            {"id": 18, "title": "greenmetrics"},
            {"id": 19, "title": "multi-greenmetrics"},
            {"id": 20, "title": "datetime"},
            {"id": 21, "title": "formatted-text-input"},
        ],
    )


def downgrade_():
    # Clean up the seeded data if we roll back
    op.execute('DELETE FROM "InputType" WHERE id BETWEEN 1 AND 21')


def upgrade_audit_log():
    pass


def downgrade_audit_log():
    pass
