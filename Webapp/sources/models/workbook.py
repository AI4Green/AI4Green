from sources.extensions import db

from .base import Model


class WorkBook(Model):
    """
    name is the name of the workbook
    group is the workgroup the workbook is associated
    users is a list of Persons that are allowed to access a workbook
    reactions is a list of reactions within the workbook
    the name and group is used to create the PrimaryKey
    """

    __tablename__ = "WorkBook"
    __table_args__ = (db.UniqueConstraint("name", "group"),)

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.Text, nullable=False)
    abbreviation = db.Column(db.Text, nullable=False)
    group = db.Column(
        db.ForeignKey("WorkGroup.id", ondelete="CASCADE"), nullable=False, index=True
    )
    time_of_creation = db.Column(db.DateTime)

    WorkGroup = db.relationship("WorkGroup", backref="WorkBook")
    reactions = db.relationship(
        "Reaction",
        backref="WorkBook",
        cascade="all, delete",
    )
    users = db.relationship("Person", secondary="Person_WorkBook")
    retrosynthesis = db.relationship("Retrosynthesis", backref="WorkBook")
