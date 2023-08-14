from marshmallow import fields
from sources.extensions import ma
from sources import models
from .UserSchema import UserSchema

class PersonSchema(ma.SQLAlchemyAutoSchema):
    """
    Person Schema
    """

    class Meta:
        model = models.Person

    user = ma.Nested(UserSchema, attribute="user")
