from sources.news_feed import news_feed_bp
from flask_login import login_required
from flask import request, flash, redirect, url_for
from sources import db
from datetime import datetime
import pytz


@news_feed_bp.route('/new_news_post', methods=['POST'])
@login_required
def new_news_post():
    # save to database
    title = request.form["newsPostTitle"]
    message = request.form["newsPostMessage"]
    current_time = datetime.now(pytz.timezone('Europe/London')).replace(tzinfo=None)
    db.NewsItem(title=title, message=message, time=current_time)
    # success and reload admin dashboard
    flash("Your message has been posted!")
    return redirect(url_for("admin_dashboard.admin_dashboard"))


@news_feed_bp.route('/delete_news_post', methods=['POST'])
@login_required
def delete_news_post():
    # get id
    item_id = request.form["post-to-delete"]
    # delete post
    db.NewsItem[item_id].delete()
    # success and reload home
    flash("This message has been deleted!")
    return redirect(url_for("main.index"))
