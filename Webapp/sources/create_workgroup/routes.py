from flask import render_template, flash, redirect, url_for, request  # renders html templates
from flask_login import login_required, current_user  # protects a view function against anonymous users
from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField, SelectField, TextAreaField
from wtforms.validators import Length
from sources.create_workgroup import create_workgroup_bp  # imports the blueprint of the dummy route
from pony.orm import select
from sources import db, auxiliary
from sources.auxiliary import get_workgroups
from datetime import datetime
from sources.email_methods import send_notification_email


class CreateWorkgroupForm(FlaskForm):
    institution = SelectField('Institution')
    workgroup = StringField('Workgroup', validators=[Length(min=4, max=320)])
    info = TextAreaField('Further Information',  render_kw={"placeholder": "I am a Principal Investigator researching...",
                                                            "rows": "4", "cols": "50"})
    submit = SubmitField('Create Workgroup', render_kw={"class": "btn btn-primary"})


@create_workgroup_bp.route('/create_workgroup', methods=['GET', 'POST'])
@login_required
def create_workgroup():
    # must be logged in
    workgroups = get_workgroups()
    form = CreateWorkgroupForm()
    # institutions removed for now
    # get all institutions (should PIs be associated with one or more institutions?)
    # institutions = list(select(u.name for u in db.Institution))
    # add to institution options
    # form.institution.choices = institutions
    # on valid submit
    if request.method == "POST":
        # institutions removed for now
        # institution = form.institution.data
        workgroup_name = auxiliary.sanitise_user_input(form.workgroup.data)
        info = form.info.data
        # find PI object from current user email
        PI = select(u for u in db.Person if u.user.email == current_user.email).first()
        # validate no other workgroups awaiting approval
        other_workgroups_waiting_approval = select(x for x in db.WorkGroup if x.approved is False and PI in x.principal_investigator)
        if other_workgroups_waiting_approval:
            flash("You cannot create another workgroup until your first workgroup has been approved")
            return render_template('create_workgroup.html', form=form, workgroups=workgroups)
        if workgroup_name == "":
            flash("Please choose a name for your Workgroup")
            return render_template('create_workgroup.html', form=form, workgroups=workgroups)
        # validates against special characters
        if not workgroup_name.replace(' ', '').replace('-', '').isalnum():
            flash('Workgroup names cannot contain special characters!')
            return render_template('create_workgroup.html', form=form, workgroups=workgroups)
        if info == "":
            flash("Please indicate why you need to set up a Workgroup")
            return render_template('create_workgroup.html', form=form, workgroups=workgroups)
        # make sure workgroup of same name does exist (add back in institution later)
        if select(u for u in db.WorkGroup if u.name.lower() == workgroup_name.lower()).first():
            flash("A Workgroup of this name already exists. Please choose a different name.")
            return render_template('create_workgroup.html', form=form, workgroups=workgroups)
        if select(u for u in db.WorkGroup_request if u.name.lower() == workgroup_name.lower()).first():
            flash("This request has already been made. You will receive a notification when your request has been considered.")
            return render_template('create_workgroup.html', form=form, workgroups=workgroups)
        # use the one default institution for now
        institution_obj = select(x for x in db.Institution if x.name == "Test User Institution").first()
        # create new workgroup and request for admin
        new_workgroup = db.WorkGroup(name=workgroup_name, institution=institution_obj, principal_investigator=PI, approved=False)
        db.WorkGroup_request(name=workgroup_name, institution=institution_obj, principal_investigator=PI, info=info,
                             time=datetime.now(), status="active", workgroup=new_workgroup)
        # send notification to admins
        admins = select(u for u in db.User if u.role.name == "Admin")[:]
        for admin in admins:
            admin_person = select(u for u in db.Person if u.user.email == admin.email).first()
            db.Notification(person=admin_person, type="New Workgroup Request",
                            info='A PI has requested a new Workgroup. Go to the <a href=/admin_dashboard '
                                 'id="admin_dashboard">admin dashboard</a> to see this request.',
                            wg_request="yes", time=datetime.now(), status="active")
            send_notification_email(admin_person)
        # flash success message and redirect to manage workgroup page
        flash("Your Workgroup has been created. It will show as under moderation until it has been approved by the site admin."
              "You will receive a notification when your request has been approved.")
        return redirect(url_for('main.index'))
    return render_template('create_workgroup.html', form=form, workgroups=workgroups)
