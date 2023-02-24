from functools import wraps
from flask import current_app, request, redirect, url_for
from flask_login import current_user, config


def login_not_allowed(func):
    """Decorator to stop logged in users accessing login/register routes"""
    @wraps(func)
    def decorated_view(*args, **kwargs):
        if request.method in config.EXEMPT_METHODS:
            return func(*args, **kwargs)
        elif current_app.config.get('LOGIN_DISABLED'):
            return func(*args, **kwargs)
        elif current_user.is_authenticated:
            return redirect(url_for('main.index'))
        return func(*args, **kwargs)
    return decorated_view
