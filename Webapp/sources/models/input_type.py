from sources.extensions import db

from .base import Model


class InputType(Model):

    """
    todo: lookup table, need population?
    """

    __tablename__ = "InputType"

    id = db.Column(db.Integer, primary_key=True)
    title = db.Column(db.String, nullable=False)

    field = db.relationship("Field", back_populates="input_type")

    def to_dict(self):
        return {
            "id": self.id,
            "title": self.title,
        }
