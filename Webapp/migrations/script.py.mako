"""${message}

Revision ID: ${up_revision}
Revises: ${down_revision | comma,n}
Create Date: ${create_date}

"""
from alembic import op
import sqlalchemy as sa

from migrations.utils import perform_migration
${imports if imports else ""}

# revision identifiers, used by Alembic.
revision = ${repr(up_revision)}
down_revision = ${repr(down_revision)}
branch_labels = ${repr(branch_labels)}
depends_on = ${repr(depends_on)}


def upgrade(engine_name: str):
    # test if the engine name is "default"
    # add target="target_name" to change the db
    # the migration is performed on
    # e.g. "audit_log"
    if not perform_migration(engine_name):
        return None
    ${upgrades if upgrades else "return None"}


def downgrade(engine_name: str):
    # test if the engine name is "default"
    # add target="target_name" to change the db
    # the migration is performed on
    # e.g. "audit_log"
    if not perform_migration(engine_name):
        return None
    ${downgrades if downgrades else "return None"}
