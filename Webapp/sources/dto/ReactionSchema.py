from typing import Any, Optional
from marshmallow import fields
from sources.extensions import ma
from sources import models
import json

from .PersonSchema import PersonSchema

class JsonField(fields.Field):
    def _serialize(self, value: str, attr, obj, **kwargs) -> Any:
        return None if value is None else json.loads(value)

    def _deserialize(self, value, attr, data, **kwargs) -> Optional[str]:
        return None if value is None else json.dumps(value)
    
  
class ReactionSchema(ma.SQLAlchemyAutoSchema):
    """
    Reaction Schema for CSV export.
    """

    class Meta:
        model = models.Reaction
        exclude = ['reaction_class', 'date_reaction', 'green_metric', 'status']

    creator = ma.Nested(PersonSchema, attribute="creator_person")
    reaction_table_data = JsonField(attribute="reaction_table_data")
    summary_table_data = JsonField(attribute="summary_table_data")
