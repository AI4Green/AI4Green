from sources.extensions import db

from .base import Model


class Section(Model):
    __tablename__ = "Section"

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String, nullable=False)

    sort_order = db.Column(db.Integer)

    section_type_id = db.Column(
        db.Integer, db.ForeignKey("SectionType.id"), nullable=False
    )
    section_type = db.relationship("SectionType", back_populates="sections")

    template_id = db.Column(db.Integer, db.ForeignKey("Template.id"))
    template = db.relationship("Template", back_populates="sections")

    fields = db.relationship("Fields", back_populates="section")
