from marshmallow import Schema, fields


class TestResponseSchema(Schema):
    message = fields.Str()
    user_status = fields.Str()
    workgroup = fields.Str()
