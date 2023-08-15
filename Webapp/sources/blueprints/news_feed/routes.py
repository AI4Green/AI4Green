from datetime import datetime

import pytz
from flask import flash, redirect, request, url_for
from flask_login import login_required

from sources import models
from sources.extensions import db

from . import news_feed_bp


@news_feed_bp.route("/new_news_post", methods=["POST"])
@login_required
def new_news_post():
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
def delete_news_post():
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
