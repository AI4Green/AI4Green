from sources.extensions import db

from .base import Model


class SectionType(Model):
    """
    todo: Maybe we dont need this?
    """

    __tablename__ = "SectionType"

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String, nullable=False)

    sections = db.relationship("Section", back_populates="section_type")

    def to_dict(self):
        return {
            "id": self.id,
            "name": self.name,
            "sections": self.sections,  # todo: is this needed?
        }
