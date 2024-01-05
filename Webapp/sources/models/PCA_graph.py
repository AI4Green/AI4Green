from sources.extensions import db

from .base import Model

# metadata = db.Model.metadata


class PCAGraph(Model):
    """
    PCA plots can be replotted from:
    r_class: reaction class selected in reaction class dropdown
    colour_selected: colour in colour dropdown
    control_points: defined points for constrained PCA

    """

    __tablename__ = "PCA"

    id = db.Column(db.Integer, primary_key=True)

    graph_name = db.Column(db.Text, nullable=False)
    description = db.Column(db.Text, nullable=False, default="")
    r_class = db.Column(db.Text, nullable=False, default="")
    colour_selected = db.Column(db.Text, nullable=False, default="CHEM21")
    descriptors = db.Column(db.Text, nullable=False, default="[]")
    control_points = db.Column(db.JSON, nullable=False, default="{}")
    graph_data = db.Column(db.JSON, nullable=False, default="{}")
    embedding_algorithm = db.Column(db.JSON, nullable=False, default="{}")

    time_of_creation = db.Column(db.DateTime, nullable=False)
    time_of_update = db.Column(db.DateTime)
    date_of_creation = db.Column(db.DateTime)

    creator = db.Column(
        db.ForeignKey("Person.id", ondelete="CASCADE"), nullable=False, index=True
    )
    creator_person = db.relationship("Person", foreign_keys=[creator])

    status = db.Column(db.Text, nullable=False, default="active")
