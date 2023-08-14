from flask_wtf import FlaskForm
from wtforms import SelectField, SubmitField


# form class to control the workgroup request
class JoinWorkgroupForm(FlaskForm):
    workgroups = SelectField("Workgroup", coerce=str)
    submit = SubmitField(
        "Request to Join Workgroup", render_kw={"class": "btn btn-primary"}
    )
