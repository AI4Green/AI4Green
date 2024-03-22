from sources.extensions import db

from .base import Model

metadata = db.Model.metadata


class Reaction(Model):
    """
    A Reaction is defined by;
    name: identifier
    date: date of creation of the reaction
    date_reaction: date when the reaction was performed
    precursor_reaction: links to the reaction in the database that preceded it
    successor_reaction: links to the reactions that follow on from this one
    green_metric: the score given for this reaction
    reactants: links to a set of reactants for this reaction
    products: links to a set of products for this reaction
    reagents: links to a set of reagents for this reaction
    solvent: links to a set of solvents for this reaction
    name_book_comp: this is the Unique PrimaryKey id for the reaction based on the
        name of the reaction and the workbook it is performed in
    """

    __tablename__ = "Reaction"
    __table_args__ = (
        db.UniqueConstraint("workbooks", "name"),
        db.UniqueConstraint("workbooks", "reaction_id"),
    )

    id = db.Column(db.Integer, primary_key=True)
    reaction_id = db.Column(db.Text, nullable=False)
    name = db.Column(db.Text, nullable=False)
    description = db.Column(db.Text, nullable=False, default="")
    reaction_class = db.Column(db.Text, nullable=False, default="")
    creator = db.Column(
        db.ForeignKey("Person.id", ondelete="CASCADE"), nullable=False, index=True
    )
    time_of_creation = db.Column(db.DateTime, nullable=False)
    time_of_update = db.Column(db.DateTime)
    date_reaction = db.Column(db.DateTime)
    green_metric = db.Column(db.Float(53))
    workbooks = db.Column(
        db.ForeignKey("WorkBook.id", ondelete="CASCADE"), nullable=False
    )
    reactants = db.Column(db.ARRAY(db.Text()), nullable=False, default="{}")  # smiles
    products = db.Column(db.ARRAY(db.Text()), nullable=False, default="{}")  # smiles
    reagents = db.Column(db.ARRAY(db.Text()), nullable=False, default="{}")  # smiles
    solvent = db.Column(
        db.ARRAY(db.Text()), nullable=False, default="{}"
    )  # primary keys
    reaction_table_data = db.Column(db.JSON, nullable=False)
    summary_table_data = db.Column(db.JSON, nullable=False)
    reaction_smiles = db.Column(db.Text, nullable=False, default="")
    complete = db.Column(db.Text, nullable=False)
    status = db.Column(db.Text, nullable=False)
    precursor_reaction = db.relationship(
        "Reaction",
        secondary="Reaction_Reaction",
        primaryjoin="Reaction.id == Reaction_Reaction.c.reaction",
        foreign_keys="[Reaction_Reaction.c.reaction]",
    )
    successor_reaction = db.relationship(
        "Reaction",
        secondary="Reaction_Reaction",
        primaryjoin="Reaction.id == Reaction_Reaction.c.reaction_2",
        foreign_keys="[Reaction_Reaction.c.reaction_2]",
        viewonly=True,
    )
    file_attachments = db.relationship("ReactionDataFile", backref="reaction_file")
    addenda = db.relationship("ReactionNote", backref="reaction_addenda")

    creator_person = db.relationship("Person", foreign_keys=[creator])
    workbook = db.relationship("WorkBook", foreign_keys=[workbooks])

    # data_export_request_id = db.Column(db.Integer, db.ForeignKey('Reaction.id'))
    data_export_request_id = db.Column(
        db.Integer, db.ForeignKey("DataExportRequest.id")
    )
    # data_export_request = db.relationship("DataExportRequest", backref="reactions")


t_Reaction_Reaction = db.Table(
    "Reaction_Reaction",
    metadata,
    db.Column(
        "reaction",
        db.ForeignKey("Reaction.id", ondelete="CASCADE"),
        primary_key=True,
        nullable=False,
    ),
    db.Column(
        "reaction_2",
        db.ForeignKey("Reaction.id", ondelete="CASCADE"),
        primary_key=True,
        nullable=False,
        index=True,
    ),
)
