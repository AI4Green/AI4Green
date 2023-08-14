from sources.extensions import db
from .base import Model


class ReactionDataFile(Model):
    __tablename__ = "ReactionDataFile"

    uuid = db.Column(db.Text, primary_key=True)
    reaction = db.Column(
        db.ForeignKey("Reaction.id", ondelete="CASCADE"), nullable=False, index=True
    )
    storage_name = db.Column(db.Text, nullable=False)
    container_name = db.Column(db.Text, nullable=False)
    display_name = db.Column(db.Text, nullable=False)
    time_of_upload = db.Column(db.DateTime, nullable=False)
    file_details = db.Column(db.JSON, nullable=False)
