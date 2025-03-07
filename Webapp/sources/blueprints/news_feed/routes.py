from datetime import datetime

import pytz
from flask import flash, redirect, request, url_for
from flask_login import login_required
from sources import models
from sources.decorators import admin_required
from sources.extensions import db

from . import news_feed_bp


@news_feed_bp.route("/new_news_post", methods=["POST"])
@login_required
@admin_required
@news_feed_bp.doc(
    security="sessionAuth",
    description="Requires login and admin user role.",
)
def new_news_post():
    """
    Create a new news post and save it to the database
    Returns:
        flask.Response: A redirect to the admin dashboard with a message indicating whether the post was successful.
    """
    # save to database
    title = request.form["newsPostTitle"]
    message = request.form["newsPostMessage"]
    current_time = datetime.now(pytz.timezone("Europe/London")).replace(tzinfo=None)
    item = models.NewsItem(title=title, message=message, time=current_time)
    db.session.add(item)
    db.session.commit()
    # success and reload admin dashboard
    flash("Your message has been posted!")
    return redirect(url_for("admin_dashboard.admin_dashboard"))


@news_feed_bp.route("/delete_news_post", methods=["POST"])
@login_required
@admin_required
@news_feed_bp.doc(
    security="sessionAuth",
    description="Requires login and admin user role.",
)
def delete_news_post():
    """
    Delete a news post from the database and feed

    Returns:
        flask.Response: A redirect to the home page with a message indicating whether the post was deleted.
    """
    # get id
    item_id = request.form["post-to-delete"]
    # delete post
    item = (
        db.session.query(models.NewsItem).filter(models.NewsItem.id == item_id).first()
    )
    db.session.delete(item)
    db.session.commit()
    # success and reload home
    flash("This message has been deleted!")
    return redirect(url_for("main.index"))
