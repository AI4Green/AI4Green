"""cheminventory_connection

Revision ID: a7c3e91d2f40
Revises: 15d154c4b3fb
Create Date: 2026-11-03 12:00:00.000000

"""
import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision = "a7c3e91d2f40"
down_revision = "15d154c4b3fb"
branch_labels = None
depends_on = None


def upgrade(engine_name):
    globals()["upgrade_%s" % engine_name]()


def downgrade(engine_name):
    globals()["downgrade_%s" % engine_name]()


def upgrade_():
    op.create_table(
        "ChemInventoryConnection",
        sa.Column("id", sa.Integer(), nullable=False),
        sa.Column("workgroup_id", sa.Integer(), nullable=False),
        sa.Column("encrypted_token", sa.Text(), nullable=False),
        sa.Column("token_hint", sa.String(length=8), nullable=False),
        sa.Column("cheminventory_user", sa.String(length=255), nullable=True),
        sa.Column("inventory_name", sa.String(length=255), nullable=True),
        sa.Column("created_by_id", sa.Integer(), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.Column("updated_at", sa.DateTime(timezone=True), nullable=False),
        sa.ForeignKeyConstraint(["workgroup_id"], ["WorkGroup.id"]),
        sa.ForeignKeyConstraint(["created_by_id"], ["User.id"]),
        sa.PrimaryKeyConstraint("id"),
        sa.UniqueConstraint("workgroup_id"),
    )


def downgrade_():
    op.drop_table("ChemInventoryConnection")


def upgrade_update():
    pass


def downgrade_update():
    pass


def upgrade_audit_log():
    pass


def downgrade_audit_log():
    pass
