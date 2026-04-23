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
            {"id": 5, "title": "ReactionScheme"},
            {"id": 6, "title": "MultiReactionScheme"},
            {"id": 7, "title": "Radio"},
            {"id": 8, "title": "Multiple"},
            {"id": 9, "title": "File"},
            {"id": 10, "title": "ImageFile"},
            {"id": 11, "title": "Header"},
            {"id": 12, "title": "ChemicalDisposalTable"},
            {"id": 13, "title": "SortableList"},
            {"id": 14, "title": "ProjectGroupPlanTable"},
            {"id": 15, "title": "ProjectGroupHazardTable"},
            {"id": 16, "title": "YieldTable"},
            {"id": 17, "title": "MultiYieldTable"},
            {"id": 18, "title": "GreenMetrics able"},
            {"id": 19, "title": "MultiGreenMetricsTable"},
            {"id": 20, "title": "DateAndTime"},
            {"id": 21, "title": "FormattedTextInput"},
        ],
    )


def downgrade_():
    # Clean up the seeded data if we roll back
    op.execute('DELETE FROM "InputType" WHERE id BETWEEN 1 AND 21')


def upgrade_audit_log():
    pass


def downgrade_audit_log():
    pass
