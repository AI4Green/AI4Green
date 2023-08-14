from flask_wtf import FlaskForm
from wtforms import PasswordField, StringField, SubmitField
from wtforms.validators import DataRequired, Email, EqualTo


# form class to control the password reset request
class ResetPasswordRequestForm(FlaskForm):
    email = StringField(
        "Email",
        validators=[DataRequired(), Email()],
        render_kw={"placeholder": "Email", "class": "form-control form-control-lg"},
    )
    submit = SubmitField(
        "Request Password Reset", render_kw={"class": "btn btn-primary"}
    )


# form class to control the password reset function
class ResetPasswordForm(FlaskForm):
    password = PasswordField(
        "Password",
        validators=[DataRequired()],
        render_kw={"placeholder": "Password", "class": "form-control form-control-lg"},
    )
    password2 = PasswordField(
        "Repeat Password",
        validators=[DataRequired(), EqualTo("password")],
        render_kw={
            "placeholder": "Repeat Password",
            "class": "form-control form-control-lg",
        },
    )
    submit = SubmitField("Reset Password", render_kw={"class": "btn btn-primary"})
