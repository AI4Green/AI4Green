from enum import Enum

from sources.extensions import db

from .base import Model


class OutcomeValueType(Enum):
    PERCENTAGE = "percentage"
    NUMERIC = "numeric"
    RATIO = "ratio"
    BOOLEAN = "boolean"
    TEXT = "text"


class ReactionResultType(Model):
    __tablename__ = "ReactionResultType"

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String, nullable=False)
    description = db.Column(db.String, nullable=False)
    outcome_type = db.Column(db.Enum(OutcomeValueType))
    default_unit = db.Column(db.String, nullable=False)

    workgroup_id = db.Column(
        db.Integer, db.ForeignKey("WorkGroup.id"), nullable=True, index=True
    )
    workgroup = db.relationship(
        "WorkGroup",
        backref="reaction_result_types",
    )

    workbook_id = db.Column(
        db.Integer, db.ForeignKey("WorkBook.id"), nullable=True, index=True
    )
    workbook = db.relationship(
        "WorkBook",
        backref="reaction_result_types",
    )

    creator_person_id = db.Column(
        db.Integer, db.ForeignKey("Person.id"), nullable=True, index=True
    )
    creator_person = db.relationship(
        "Person",
        backref="reaction_result_types",
    )
