from sources.extensions import db

from .base import Model


class Section(Model):
    __tablename__ = "Section"

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String, nullable=False)

    sort_order = db.Column(db.Integer)

    template_id = db.Column(db.Integer, db.ForeignKey("Template.id"))
    template = db.relationship("Template", back_populates="sections")

    fields = db.relationship("Field", back_populates="section")
