from sources.extensions import ma
from sources import models

from .PersonSchema import PersonSchema

class ReactionNoteSchema(ma.SQLAlchemyAutoSchema):
    """
    ReactionNote Schema for serialisation.
    """

    class Meta:
        model = models.ReactionNote
        datetimeformat = '%Y-%m-%d %H:%M:%S'

    author = ma.Nested(PersonSchema, attribute="Author")
