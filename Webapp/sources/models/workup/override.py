from sources.extensions import db
from sources.models.base import Model


class InstanceSectionOverride(Model):
    """
    How one instance shows one of its sections (from the template or added by the
    instance): its position and whether it is hidden.

    Rows are only written once a user reorders or hides steps in an instance, so an
    unchanged instance has none. Template sections and fields are never copied.
    """

    __tablename__ = "InstanceSectionOverride"

    template_instance_id = db.Column(
        db.Integer,
        db.ForeignKey("WorkupInstance.id", ondelete="CASCADE"),
        primary_key=True,
    )
    section_id = db.Column(
        db.Integer,
        db.ForeignKey("Section.id", ondelete="CASCADE"),
        primary_key=True,
    )
    sort_order = db.Column(db.Integer, nullable=False)
    is_hidden = db.Column(db.Boolean, nullable=False, default=False)

    template_instance = db.relationship(
        "WorkupInstance", back_populates="section_overrides"
    )
    section = db.relationship("Section")
