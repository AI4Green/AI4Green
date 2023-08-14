from datetime import datetime

from sources.extensions import db

from .base import Model

metadata = db.Model.metadata


class Person(Model):
    """
    A Person is a user, and can have some attributes that give them further access
    rights within the Reaction Database

    email: Required unique identifier of the user
    fullname: Required
    workgroups: Optional Set of links to the workgroups (labs) that this user manages
    workgroup_principal_investigator: This is a set of links to the workgroups which this user is the PI/owner of
    workgroup_senior_researcher: This is a set of links to the workgroups which this user is a senior researcher in
    workgroup_standard_member: This is a set of links to the workgroups which this user is a standard member in
    workbook_user: This is a Set of links to the workbooks this user has access to
    """

    __tablename__ = "Person"

    def __init__(self) -> None:
        self.time_of_creation = datetime.now()
        super().__init__()

    id = db.Column(db.Integer, primary_key=True)
    time_of_creation = db.Column(db.DateTime)

    user = db.relationship("User", backref="Person", uselist=False)
    workgroup_principal_investigator = db.relationship(
        "WorkGroup",
        secondary="Person_WorkGroup",
        primaryjoin="Person.id == Person_WorkGroup.c.person",
        secondaryjoin="WorkGroup.id == Person_WorkGroup.c.workgroup",
        lazy="dynamic",
    )  # PI is also the owner
    workgroup_senior_researcher = db.relationship(
        "WorkGroup",
        secondary="Person_WorkGroup_2",
        primaryjoin="Person.id == Person_WorkGroup_2.c.person",
        secondaryjoin="WorkGroup.id == Person_WorkGroup_2.c.workgroup",
        lazy="dynamic",
    )
    workgroup_standard_member = db.relationship(
        "WorkGroup",
        secondary="Person_WorkGroup_3",
        primaryjoin="Person.id == Person_WorkGroup_3.c.person",
        secondaryjoin="WorkGroup.id == Person_WorkGroup_3.c.workgroup",
        lazy="dynamic",
    )

    workgroup_principal_investigator_request = db.relationship(
        "WorkGroupRequest"
    )  # PI is also the owner
    person_notification = db.relationship("Notification")
    person_status_change = db.relationship(
        "WGStatusRequest", foreign_keys="[WGStatusRequest.person]"
    )
    person_status_change_pi = db.relationship(
        "WGStatusRequest", foreign_keys="[WGStatusRequest.principal_investigator]"
    )
    person_status_change_wb = db.relationship(
        "WBStatusRequest", foreign_keys="[WBStatusRequest.person]"
    )
    person_status_change_wb_sr_pi = db.relationship(
        "WBStatusRequest", foreign_keys="[WBStatusRequest.pi_sr]"
    )
    workbook_user = db.relationship("WorkBook", secondary="Person_WorkBook")
    reaction = db.relationship("Reaction")
    reactionNote = db.relationship("ReactionNote")


t_Person_WorkGroup = db.Table(
    "Person_WorkGroup",
    metadata,
    db.Column(
        "person",
        db.ForeignKey("Person.id", ondelete="CASCADE"),
        primary_key=True,
        nullable=False,
    ),
    db.Column(
        "workgroup",
        db.ForeignKey("WorkGroup.id", ondelete="CASCADE"),
        primary_key=True,
        nullable=False,
        index=True,
    ),
)


t_Person_WorkGroup_2 = db.Table(
    "Person_WorkGroup_2",
    metadata,
    db.Column(
        "person",
        db.ForeignKey("Person.id", ondelete="CASCADE"),
        primary_key=True,
        nullable=False,
    ),
    db.Column(
        "workgroup",
        db.ForeignKey("WorkGroup.id", ondelete="CASCADE"),
        primary_key=True,
        nullable=False,
        index=True,
    ),
)


t_Person_WorkGroup_3 = db.Table(
    "Person_WorkGroup_3",
    metadata,
    db.Column(
        "person",
        db.ForeignKey("Person.id", ondelete="CASCADE"),
        primary_key=True,
        nullable=False,
    ),
    db.Column(
        "workgroup",
        db.ForeignKey("WorkGroup.id", ondelete="CASCADE"),
        primary_key=True,
        nullable=False,
        index=True,
    ),
)


t_Person_WorkBook = db.Table(
    "Person_WorkBook",
    metadata,
    db.Column(
        "person",
        db.ForeignKey("Person.id", ondelete="CASCADE"),
        primary_key=True,
        nullable=False,
    ),
    db.Column(
        "workbook",
        db.ForeignKey("WorkBook.id", ondelete="CASCADE"),
        primary_key=True,
        nullable=False,
        index=True,
    ),
)
