import uuid
from datetime import datetime, timedelta
from enum import Enum

import pytz
from sources.extensions import db

from .base import Model


class ApprovalStatus(Enum):
    PENDING = "PENDING"
    APPROVED = "APPROVED"
    REJECTED = "REJECTED"
    EXPIRED = "EXPIRED"


# Association table for the many-to-many relationship between Person and DataExportRequests
reaction_approval_request_approvers = db.Table(
    "reaction_approval_request_approvers",
    db.Model.metadata,
    db.Column(
        "reaction_approval_request_id",
        db.Integer,
        db.ForeignKey("DataExportRequest.id"),
    ),
    db.Column("person_id", db.Integer, db.ForeignKey("Person.id")),
    db.Column("approved", db.Boolean, default=False),
    db.Column("responded", db.Boolean, default=False),
)


class ReactionApprovalRequest(Model):
    """
    A Data export request is made for a workbook. It will export the contents in a user's chosen format
    """

    __tablename__ = "DataExportRequest"

    def __init__(self, **kwargs) -> None:
        self.time_of_request = datetime.now(pytz.timezone("Europe/London")).replace(
            tzinfo=None
        )
        self.uuid = str(uuid.uuid4())
        super().__init__(**kwargs)

    time_of_request = db.Column(
        db.DateTime,
        nullable=False,
    )

    id = db.Column(db.Integer, primary_key=True)

    requestor = db.Column(
        db.ForeignKey("Person.id", ondelete="CASCADE"), nullable=False, index=True
    )
    requestor_person = db.relationship("Person", foreign_keys=[requestor])

    reaction = db.Column(
        db.ForeignKey("Reaction.id", ondelete="CASCADE"), nullable=False, index=True
    )
    Reaction = db.relationship(
        "Reaction",
        backref="reaction_approval_request",
    )

    workbook = db.Column(
        db.ForeignKey("WorkBook.id", ondelete="CASCADE"), nullable=False, index=True
    )
    WorkBook = db.relationship(
        "WorkBook",
        backref="reaction_approval_requests",
    )

    workgroup = db.Column(
        db.ForeignKey("WorkGroup.id", ondelete="CASCADE"), nullable=False, index=True
    )
    WorkGroup = db.relationship(
        "WorkGroup", backref="reaction_approval_request", cascade="all, delete"
    )
    institution = db.Column(
        db.ForeignKey("Institution.id", ondelete="CASCADE"), nullable=False, index=True
    )
    Institution = db.relationship("Institution", backref="reaction_approval_request")

    # Establish a many-to-many relationship with Person for required_approvers
    required_approvers = db.relationship(
        "Person",
        secondary=reaction_approval_request_approvers,
        backref="reaction_approval_requests",
    )

    status = db.Column(db.Enum(ApprovalStatus), default=ApprovalStatus.PENDING.value)
