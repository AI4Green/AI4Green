from marshmallow import fields, INCLUDE
from sources.extensions import ma
from sources import models

class UserSchema(ma.SQLAlchemyAutoSchema):
    """
    User Schema
    """

    class Meta:
        model = models.User
        fields = ['fullname']

