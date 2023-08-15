from flask_wtf import FlaskForm
from wtforms import PasswordField, StringField, SubmitField
from wtforms.validators import DataRequired, Email, EqualTo


# form class to ask for new password
class UpdateEmailForm(FlaskForm):
    old_password = PasswordField(
        "Password",
        validators=[DataRequired()],
        render_kw={"placeholder": "Password", "class": "form-control form-control-lg"},
    )
    email = StringField(
        "New Email",
        validators=[DataRequired(), Email()],
        render_kw={"placeholder": "New Email", "class": "form-control form-control-lg"},
    )
    email2 = StringField(
        "Repeat New Email",
        validators=[DataRequired(), EqualTo("email"), Email()],
        render_kw={
            "placeholder": "Repeat New Email",
            "class": "form-control form-control-lg",
        },
    )
    submit = SubmitField("Update Email", render_kw={"class": "btn btn-primary"})
