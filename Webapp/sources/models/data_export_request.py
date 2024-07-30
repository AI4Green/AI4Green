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


class ExportFormat(Enum):
    RDF = "RDF"
    PDF = "PDF"
    ELN = "ELN"
    SURF = "SURF"
    CSV = "CSV"
    JSON = "JSON"
    SI = "SI"


# Association table for the many-to-many relationship between Person and DataExportRequests
data_export_request_approvers = db.Table(
    "data_export_request_approvers",
    db.Model.metadata,
    db.Column(
        "data_export_request_id", db.Integer, db.ForeignKey("DataExportRequest.id")
    ),
    db.Column("person_id", db.Integer, db.ForeignKey("Person.id")),
    db.Column("approved", db.Boolean, default=False),
    db.Column("responded", db.Boolean, default=False),
)

# Association table for the many-to-many relationship between DataExportRequest and WorkBook
data_export_request_workbooks = db.Table(
    "data_export_request_workbook_association",
    db.Model.metadata,
    db.Column(
        "data_export_request_id", db.Integer, db.ForeignKey("DataExportRequest.id")
    ),
    db.Column("workbook_id", db.Integer, db.ForeignKey("WorkBook.id")),
)

# Association table for the many-to-many relationship between DataExportRequest and Reaction
data_export_request_reactions = db.Table(
    "data_export_request_reactions",
    db.Model.metadata,
    db.Column("data_export_request_id", db.Integer, db.ForeignKey("DataExportRequest.id")),
    db.Column("reaction_id", db.Integer, db.ForeignKey("Reaction.id"))
)


class DataExportRequest(Model):
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
    data_format = db.Column(db.Enum(ExportFormat))

    requestor = db.Column(
        db.ForeignKey("Person.id", ondelete="CASCADE"), nullable=False, index=True
    )
    requestor_person = db.relationship("Person", foreign_keys=[requestor])

    # Establish a many-to-many relationship with Person for required_approvers
    required_approvers = db.relationship(
        "Person",
        secondary=data_export_request_approvers,
        backref="data_export_requests",
    )

    status = db.Column(db.Enum(ApprovalStatus), default=ApprovalStatus.PENDING.value)

    uuid = db.Column(db.Text)  # Unique container name
    hash = db.Column(db.Text)  # to confirm zip download contents

    reactions = db.relationship(
        "Reaction",
        secondary=data_export_request_reactions,
        backref="data_export_request"
    )

    # supports multiple workbooks if desired in future functionality
    workbooks = db.relationship(
        "WorkBook",
        secondary=data_export_request_workbooks,
        backref="data_export_requests",
    )

    workgroup = db.Column(
        db.ForeignKey("WorkGroup.id", ondelete="CASCADE"), nullable=False, index=True
    )
    WorkGroup = db.relationship(
        "WorkGroup", backref="data_export_request", cascade="all, delete"
    )
    institution = db.Column(
        db.ForeignKey("Institution.id", ondelete="CASCADE"), nullable=False, index=True
    )
    Institution = db.relationship("Institution", backref="data_export_request")
