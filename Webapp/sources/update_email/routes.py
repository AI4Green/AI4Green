from flask import render_template, flash, redirect, url_for, request  # renders html templates
from flask_login import login_required  # protects a view function against anonymous users
from sources.email_methods import send_password_reset_email
from sources.update_email import update_email_bp  # imports the blueprint of the dummy route
from sources.update_email.forms import UpdateEmailForm
from flask_login import current_user
from pony.orm import select
from sources import db
from sources.auxiliary import get_workgroups, get_notification_number


# Go to the update email/password page
@update_email_bp.route('/update_email_password', methods=['GET', 'POST'])
@login_required
def update_email_password():
    # must be logged in
    current_email = current_user.email
    user = select(u for u in db.User if current_user.username == u.username).first()
    if request.method == 'POST':
        if request.form.get('ChangeEmail') == 'Value1':
            return redirect(url_for('update_email.update_email'))
        elif request.form.get('UpdatePassword') == 'Value2':
            send_password_reset_email(user)
            flash('Please check your email for instructions on how to change your password')
            return redirect(url_for('main.index'))
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    return render_template('update_email_password.html', email=current_email, workgroups=workgroups,
                           notification_number=notification_number)


# Go to the update email page
@update_email_bp.route('/update_email', methods=['GET', 'POST'])
@login_required
def update_email():
    # must be logged in
    form = UpdateEmailForm()
    if form.validate_on_submit():
        # if authentication is successful
        user = select(u for u in db.User if current_user.username == u.username).first()
        if user is None or not user.check_password(form.old_password.data):
            flash('Invalid password')
            return redirect(url_for('update_email.update_email'))
        # check email is not already registered
        if form.email.data in list(select(u.email for u in db.User)):
            flash('This email is already registered to an account')
            return redirect(url_for('update_email.update_email'))
        # change user's email
        user.email = form.email.data
        # alert user of success and redirect to home
        flash("Your email has been successfully updated!")
        return redirect(url_for('main.index'))
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    return render_template('update_email.html', form=form, workgroups=workgroups,
                           notification_number=notification_number)

