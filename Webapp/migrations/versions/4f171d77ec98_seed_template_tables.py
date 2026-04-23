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
            {"id": 1, "title": "text"},  # single line text
            {"id": 2, "title": "textarea"},  # long text
            {"id": 3, "title": "number"},  # numeric values
            {"id": 4, "title": "select"},  # single choice dropdown
            {"id": 5, "title": "multiselect"},  # multiple choice
            {"id": 6, "title": "checkbox"},  # boolean true/false
            {"id": 7, "title": "date"},  # dates (review, assessment)
            {"id": 8, "title": "file"},  # SDS uploads, documents
            {"id": 9, "title": "url"},  # references / external links
            {"id": 10, "title": "richtext"},  # formatted procedures
            {"id": 11, "title": "repeatable_group"},  # repeated structured entries
        ],
    )


def downgrade_():
    # Clean up the seeded data if we roll back
    op.execute('DELETE FROM "InputType" WHERE id BETWEEN 1 AND 11')


def upgrade_audit_log():
    pass


def downgrade_audit_log():
    pass
