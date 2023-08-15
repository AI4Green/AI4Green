from functools import wraps

from flask import current_app, redirect, request, url_for
from flask_login import config, current_user

from sources import login, models


def login_not_allowed(func):
    """Decorator to stop logged in users accessing login/register routes"""

    @wraps(func)
    def decorated_view(*args, **kwargs):
        if request.method in config.EXEMPT_METHODS:
            return func(*args, **kwargs)
        elif current_app.config.get("LOGIN_DISABLED"):
            return func(*args, **kwargs)
        elif current_user.is_authenticated:
            return redirect(url_for("main.index"))
        return func(*args, **kwargs)

    return decorated_view


@login.user_loader
def load_user(user_id) -> models.User:
    """Flask-Login keeps track of the logged in user
    by storing its unique identifier in Flask's user
    session. Each time the logged-in user navigates
    to a new page, Flask-Login retrieves the ID of
    the user from the session, and then loads that
    user into memory by the user loader function"""
    return models.User.query.get(user_id)
