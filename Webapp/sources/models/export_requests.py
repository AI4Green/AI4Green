from sources.extensions import db

from .base import Model


class DataExportRequest(Model):
    """
    A Data export request is made for a workbook. It will export the contents in a user's chosen format

    """

    __tablename__ = "WorkGroup_request"
    __table_args__ = (db.UniqueConstraint("name", "institution"),)

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.Text, nullable=False)
    principal_investigator = db.Column(
        db.ForeignKey("Person.id", ondelete="CASCADE"), nullable=False, index=True
    )
    info = db.Column(db.Text, nullable=False)
    time = db.Column(db.DateTime, nullable=False)
    status = db.Column(db.Text, nullable=False)
    institution = db.Column(
        db.ForeignKey("Institution.id", ondelete="CASCADE"), nullable=False, index=True
    )
    workgroup = db.Column(
        db.ForeignKey("WorkGroup.id", ondelete="CASCADE"), nullable=False, index=True
    )

    pi = db.relationship(
        "Person", foreign_keys=[principal_investigator], backref="WorkGroupRequest"
    )
    Institution = db.relationship("Institution", backref="WorkGroupRequest")
    WorkGroup = db.relationship(
        "WorkGroup", backref="WorkGroupRequest", cascade="all, delete"
    )
