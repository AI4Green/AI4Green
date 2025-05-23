"""
This module stores web form classes for user authentication
"""
from flask_wtf import FlaskForm  # imports the FlaskForm base class
from flask_wtf.recaptcha import RecaptchaField
from sources import models

# The validators argument is used to attach validation behaviors to fields.
# The ValidationError validator triggers a validation error if user data
# are already in the database.
# The DataRequired validator checks that the field is not submitted empty.
# The Email validator ensures that what the user types in the email field.
# matches the structure of an email address.
# The EqualTo validator makes sure that the value for the second password field
# is identical to the one for the first password field.
from sources.extensions import db  # imports the User class from sources/models.py
from sqlalchemy import func

# The Flask-WTF extension uses Python classes to represent web forms
# A form class defines the fields of the form as class variables
from wtforms import BooleanField, EmailField, PasswordField, StringField, SubmitField

# imports the classes representing the field types from the WTForms package,
# since the Flask-WTF extension does not provide customized versions
from wtforms.validators import (
    DataRequired,
    Email,
    EqualTo,
    Length,
    Regexp,
    ValidationError,
)


class LoginForm(FlaskForm):  # this class defines the login form fields
    username = StringField(
        "Username or Email",
        validators=[DataRequired("Username or email is required")],
        render_kw={
            "placeholder": "Username or email",
            "class": "form-control form-control-lg",
        },
    )
    password = PasswordField(
        "Password",
        validators=[DataRequired("Password is required")],
        render_kw={"placeholder": "Password", "class": "form-control form-control-lg"},
    )
    remember_me = BooleanField("Remember Me")
    submit = SubmitField("Sign In", render_kw={"class": "btn btn-primary"})


class RegistrationForm(FlaskForm):  # this class defines the registration form fields
    username = StringField(
        "Username",
        validators=[
            DataRequired("Username is required"),
            Length(min=4, max=32),
            Regexp("^[a-zA-Z0-9_]+$", message="Invalid character included in username"),
        ],
        render_kw={"placeholder": "Username", "class": "form-control form-control-lg"},
    )
    email = EmailField(
        "Email",
        validators=[DataRequired("Email is required"), Email(), Length(max=320)],
        render_kw={"placeholder": "Email", "class": "form-control form-control-lg"},
    )
    fullname = StringField(
        "Full Name",
        validators=[DataRequired("Full name is required"), Length(max=50)],
        render_kw={"placeholder": "Full Name", "class": "form-control form-control-lg"},
    )
    password = PasswordField(
        "Password",
        validators=[DataRequired("Password is required."), Length(min=8)],
        render_kw={"placeholder": "Password", "class": "form-control form-control-lg"},
    )
    password2 = PasswordField(
        "Confirm Password",
        # validators=[DataRequired(), EqualTo("password")],
        validators=[
            DataRequired(message="Confirm Password is required."),
            EqualTo(
                "password",
                message="The passwords do not match. Please ensure both password fields are identical.",
            ),
        ],
        render_kw={
            "placeholder": "Repeat Password",
            "class": "form-control form-control-lg",
        },
    )

    # a user is asked to type the password twice to reduce the risk of a typo
    privacy = BooleanField(
        "I have read and agreed to the conditions in the privacy notice",
        render_kw={"onclick": "enableSubmit()"},
    )
    hazard_disclaimer = BooleanField(
        "I have read and understood the hazard disclaimer",
        render_kw={"onclick": "enableSubmit()"},
    )

    recaptcha = RecaptchaField()

    submit = SubmitField(
        "Register", render_kw={"class": "btn btn-primary", "disabled": "true"}
    )

    """WTForms takes methods like validate_<field_name> as custom validators
    and invokes them in addition to the stock validators. The validators
    validate_username and validate_email are used to make sure that the
    username and email address entered by the user are not already in
    the database, so these two methods issue database queries expecting
    there will be no results. In the event a result exists, a validation
    error is triggered by raising ValidationError. The message included
    as the argument in the exception will be the message that will be
    displayed next to the field for the user to see."""

    def validate_username(self, username: StringField) -> None:
        user = (
            db.session.query(models.User)
            .filter(func.lower(models.User.username) == username.data.lower())
            .first()
        )
        if user is not None:
            raise ValidationError(
                "This username is taken. Please choose a different username."
            )

    def validate_email(self, email: StringField) -> None:
        user = (
            db.session.query(models.User)
            .filter(func.lower(models.User.email) == email.data.lower())
            .first()
        )
        if user is not None:
            raise ValidationError(
                "A user is already registered with this email address. Try logging in or resetting your password."
            )

    def validate_privacy(self, privacy: BooleanField) -> None:
        if privacy.data is not True:
            raise ValidationError(
                "You need to agree to the privacy notice in order to registration"
            )

    def validate_hazard_disclaimer(self, hazard_disclaimer: BooleanField) -> None:
        if hazard_disclaimer.data is not True:
            raise ValidationError(
                "You need to agree to the hazard disclaimer in order to registration"
            )
