from sources.extensions import db

from .base import Model


class Section(Model):
    __tablename__ = "Section"

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String, nullable=False)

    sort_order = db.Column(db.Integer)

    template_id = db.Column(db.Integer, db.ForeignKey("COSHHTemplate.id"))
    template = db.relationship("COSHHTemplate", back_populates="sections")

    fields = db.relationship("Field", back_populates="section")

    section_type = db.Column(db.String(50), nullable=False)

    __mapper_args__ = {
        "polymorphic_on": section_type,
        "polymorphic_identity": "section",
    }

    def to_dict(self):
        return {
            "id": self.id,
            "name": self.name,
            "sortOrder": self.sort_order,
            "templateId": self.template_id,
            "fields": [x.to_dict() for x in self.fields],
        }
