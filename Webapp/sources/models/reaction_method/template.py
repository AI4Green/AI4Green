from typing import Any, Dict

from sources.extensions import db

from ..base import Model
from .method import _now


class ReactionMethodTemplate(Model):
    """
    A saved reaction method, shared within one workbook.

    `setup_conditions` holds the method's setup and conditions as serialised by
    ReactionMethod.setup_dict and conditions_dict. `additions` holds the order of
    addition as a pattern: components are recorded by role and position in the
    reaction table (e.g. the second reagent), not by SMILES, so the template can be
    applied to reactions with different compounds.
    """

    __tablename__ = "ReactionMethodTemplate"

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.Text, nullable=False)
    description = db.Column(db.Text, nullable=False, default="")
    workbook_id = db.Column(
        db.Integer,
        db.ForeignKey("WorkBook.id", ondelete="CASCADE"),
        nullable=False,
        index=True,
    )
    workbook = db.relationship(
        "WorkBook",
        backref=db.backref(
            "reaction_method_templates",
            cascade="all, delete-orphan",
            passive_deletes=True,
        ),
    )
    creator_id = db.Column(
        db.Integer, db.ForeignKey("User.id", ondelete="SET NULL"), nullable=True
    )
    creator = db.relationship("User")
    time_of_creation = db.Column(db.DateTime, nullable=False, default=_now)
    setup_conditions = db.Column(db.JSON, nullable=False, default=dict)
    additions = db.Column(db.JSON, nullable=False, default=list)

    def to_dict(self) -> Dict[str, Any]:
        """
        Serialises the template.

        Returns:
            Dictionary representation of the template.
        """
        return {
            "id": self.id,
            "name": self.name,
            "description": self.description,
            "creatorId": self.creator_id,
            "creator": self.creator.fullname if self.creator else None,
            "setupConditions": self.setup_conditions or {},
            "additions": self.additions or [],
        }
