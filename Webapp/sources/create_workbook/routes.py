from flask import render_template, flash, redirect, url_for, request  # renders html templates
from flask_login import login_required, current_user  # protects a view function against anonymous users
from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField
from wtforms.validators import Length
from sources.create_workbook import create_workbook_bp  # imports the blueprint of the dummy route
from pony.orm import select
from sources import db, auxiliary
from sources.auxiliary import get_workgroups, get_notification_number, security_pi_sr_workgroup


class CreateWorkbookForm(FlaskForm):
    """This form is used for creating a workbook"""
    workbook = StringField('Workbook', validators=[Length(min=4, max=320)])
    abbreviation = StringField('3 letter abbreviation', validators=[Length(min=3, max=3)])
    submit = SubmitField('Create Workbook', render_kw={"class": "btn btn-primary"})


@create_workbook_bp.route('/create_workbook/<workgroup>', methods=['GET', 'POST'])
@login_required
def create_workbook(workgroup):
    # must be logged in and a SR or PI of the workgroup
    if not security_pi_sr_workgroup(workgroup):
        flash("You do not have permission to view this page")
        return redirect(url_for('main.index'))
    """Creates a workbook using a FlaskForm. The book is created by a PI or SR and they must provide the name
    and workgroup the workbook should belong to"""
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    form = CreateWorkbookForm()
    wg = select(x for x in db.WorkGroup if x.name == workgroup).first()
    # on valid submit
    if request.method == "POST" and form.validate_on_submit():
        # validates against special characters in workbook name
        workbook_name = auxiliary.sanitise_user_input(form.workbook.data)
        if not workbook_name.replace(' ', '').replace('-', '').isalnum():
            flash('Workgroup names cannot contain special characters!')
            return render_template('create_workgroup.html', form=form, workgroups=workgroups)
        # make sure workbook of same name does exist within workbook
        if select(u for u in db.WorkBook if u.group == wg and u.name.lower() == workbook_name.lower()).first():
            flash("A workbook of this name already exists!")
            return render_template('create_workbook.html', form=form, workgroups=workgroups,
                                   notification_number=notification_number, workgroup=workgroup)
        # validate workbook abbreviation
        workbook_abbreviation = auxiliary.sanitise_user_input(form.abbreviation.data)
        if not workbook_abbreviation.isalnum():
            flash("Workbook abbreviation must only consist of alphanumeric characters!")
            return render_template('create_workbook.html', form=form, workgroups=workgroups,
                                   notification_number=notification_number, workgroup=workgroup)
        if select(u for u in db.WorkBook if u.group == wg and u.abbreviation.lower() == workbook_abbreviation.lower()).first():
            flash("A workbook with this abbreviation already exists!")
            return render_template('create_workbook.html', form=form, workgroups=workgroups,
                                   notification_number=notification_number, workgroup=workgroup)
        # find workgroup and user object
        workgroup_obj = db.WorkGroup[wg.id]
        user_obj = select(u for u in db.Person if u.user.email == current_user.email).first()
        # Create new workbook
        db.WorkBook(name=workbook_name, abbreviation=workbook_abbreviation, group=workgroup_obj, users=user_obj)
        wb_obj = select(x for x in db.WorkBook if x.group.name == workgroup and x.name == workbook_name).first()
        # flash success message and redirect to manage workgroup page
        flash("Workbook has been created!")
        return redirect(url_for('manage_workbook.manage_workbook', workgroup=workgroup, has_request="no",
                                workbook=wb_obj.id))
    else:
        return render_template('create_workbook.html', form=form, workgroups=workgroups,
                           notification_number=notification_number, workgroup=workgroup)
