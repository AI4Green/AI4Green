from sources.extensions import db

from .base import Model


class Comment(Model):
    __tablename__ = "Comment"

    id = db.Column(db.Integer, primary_key=True)
    value = db.Column(db.Text, nullable=False)
    time_of_comment = db.Column(db.DateTime, nullable=False)

    user_id = db.Column(db.Integer, db.ForeignKey("User.id"), nullable=False)
    user = db.relationship("User", backref="comment")

    field_response_id = db.Column(
        db.Integer, db.ForeignKey("FieldResponse.id"), nullable=False
    )
    field_response = db.relationship("FieldResponse", back_populates="comment")

    read = db.Column(db.Boolean, nullable=False, default=False)
