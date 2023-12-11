from datetime import datetime

import pytz
from flask_login import UserMixin
from sources.extensions import db
from werkzeug.security import check_password_hash, generate_password_hash

from .base import Model


class User(Model, UserMixin):
    """The User class inherits from db.entity which is used as a base class
    for all entities stored in the database and UserMixin which provides
    default implementations for all of these properties and methods
    This class defines attributes as class variables where the kind
    and type of attribute and attribute options are defined."""

    __tablename__ = "User"

    def __init__(self, **kwargs) -> None:
        self.time_of_creation = datetime.now(pytz.timezone("Europe/London")).replace(
            tzinfo=None
        )
        super().__init__(**kwargs)

    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.Text, nullable=False, unique=True)
    email = db.Column(db.Text, nullable=False, unique=True)
    fullname = db.Column(db.Text, nullable=False)
    password_hash = db.Column(db.Text, nullable=False)
    role = db.Column(
        db.ForeignKey("Role.id", ondelete="CASCADE"),
        nullable=False,
        index=True,
        default=2,
    )
    person = db.Column(db.ForeignKey("Person.id"), nullable=False, index=True)
    hazard_colors = db.Column(
        db.JSON,
        default={
            "Recommended": "#00ff00",
            "Problematic": "#ffff00",
            "Hazardous": "#ff0000",
            "HighlyHazardous": "#8B0000",
            "Recommended_text": "#000000",
            "Problematic_text": "#000000",
            "Hazardous_text": "#000000",
            "HighlyHazardous_text": "#ffffff",
        },
        nullable=False,
    )
    time_of_creation = db.Column(db.DateTime)
    Role = db.relationship("Role")

    """Password hashing is implemented by the two following methods"""

    def check_password(self, password):
        return check_password_hash(self.password_hash, password)

    @classmethod
    def set_password(cls, password):
        return generate_password_hash(password)
