from sources.extensions import db

from .base import Model


class WorkGroup(Model):
    """
    A workgroup is a Person's lab - if they run one.
    A workgroup thus has a name and person associated with it, that gives the Primary key
    A workgroup can have other owners - this makes sense for larger organisations
    workgroup_members is the list of Person's who are part of this workgroup
    workbooks are a set of projects (lab books) that the lab maintains
    """

    __tablename__ = "WorkGroup"
    __table_args__ = (db.UniqueConstraint("name", "institution"),)

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.Text, nullable=False)
    institution = db.Column(
        db.ForeignKey("Institution.id", ondelete="CASCADE"), nullable=False, index=True
    )
    approved = db.Column(db.Boolean)
    time_of_creation = db.Column(db.DateTime)

    WorkGroup_Status_Request = db.relationship("WGStatusRequest")
    workbooks = db.relationship("WorkBook", cascade="all, delete")
    request = db.relationship("WorkGroupRequest", cascade="all, delete")
    principal_investigator = db.relationship(
        "Person",
        secondary="Person_WorkGroup",
        backref=db.backref("principal_workgroups", lazy="dynamic"),
    )
    senior_researcher = db.relationship(
        "Person",
        secondary="Person_WorkGroup_2",
        backref=db.backref("senior_workgroups", lazy="dynamic"),
    )
    standard_member = db.relationship(
        "Person",
        secondary="Person_WorkGroup_3",
        backref=db.backref("standard_workgroups", lazy="dynamic"),
    )
    Institution = db.relationship("Institution", backref="WorkGroup")
