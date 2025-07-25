"""email_verification

Revision ID: c57d03c696a9
Revises: 436a03e01053
Create Date: 2024-04-22 13:49:05.529796

"""

from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = "c57d03c696a9"
down_revision = "436a03e01053"
branch_labels = None
depends_on = None


def upgrade(engine_name: str):
    if engine_name == "audit_log":
        return None
    # ### commands auto generated by Alembic - please adjust! ###
    with op.batch_alter_table("User", schema=None) as batch_op:
        batch_op.add_column(
            sa.Column(
                "is_verified",
                sa.Boolean(),
                server_default=sa.text("false"),
                nullable=True,
            )
        )
        batch_op.add_column(sa.Column("verified_on", sa.DateTime(), nullable=True))

    # ### end Alembic commands ###


def downgrade(engine_name: str):
    if engine_name == "audit_log":
        return None
    # ### commands auto generated by Alembic - please adjust! ###
    with op.batch_alter_table("User", schema=None) as batch_op:
        batch_op.drop_column("verified_on")
        batch_op.drop_column("is_verified")

    # ### end Alembic commands ###
