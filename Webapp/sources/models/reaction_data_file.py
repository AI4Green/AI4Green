import sqlalchemy as sa
from sources.extensions import db

from .base import Model


class ReactionDataFile(Model):
    __tablename__ = "ReactionDataFile"

    uuid = db.Column(db.Text, primary_key=True)
    reaction = db.Column(
        db.ForeignKey("Reaction.id", ondelete="CASCADE"), nullable=False, index=True
    )
    storage_name = db.Column(db.Text, nullable=False)  # Azure storage property
    container_name = db.Column(db.Text, nullable=False)  # Azure storage property
    display_name = db.Column(db.Text, nullable=False)  # Name the user sees
    time_of_upload = db.Column(db.DateTime, nullable=False)
    file_details = db.Column(
        db.JSON, nullable=False
    )  # size, mimetype, and file extension
    sha256_checksum = db.Column(db.Text, nullable=False, server_default="not_available")
    autogenerated = db.Column(
        db.Boolean, nullable=False, server_default=sa.sql.false(), default=False
    )
