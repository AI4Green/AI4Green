from flask import render_template
from sources import app, db
import jwt
from time import time
from flask_mail import Message
from sources import mail


# method to send an email given all the information
def send_email(subject, sender, recipients, text_body, html_body):
    msg = Message(subject, sender=sender, recipients=recipients)
    msg.body = text_body
    msg.html = html_body
    mail.send(msg)


# get token with expiry time
def get_reset_password_token(user):
    expires_in = 600
    return jwt.encode(
        {'reset_password': user.id, 'exp': time() + expires_in},
        app.config['SECRET_KEY'], algorithm='HS256')


# send reset email
def send_password_reset_email(user):
    token = get_reset_password_token(user)
    send_email('AI4Chem Reset Your Password',
               sender="admin@ai4green.app",
               recipients=[user.email],
               text_body=render_template('email/reset_password_text.txt',
                                         user=user, token=token),
               html_body=render_template('email/reset_password_text.html',
                                         user=user, token=token))


# send reset email
def send_notification_email(person):
    send_email('You have a new AI4Green notification',
               sender="admin@ai4green.app",
               recipients=[person.user.email],
               text_body=render_template('email/notification_text.txt', user=person.user),
               html_body=render_template('email/notification_text.html', user=person.user))


# send reset email
def send_password_reset_email_test(user):
    token = get_reset_password_token(user)
    return render_template('email/reset_password_text.html', user=user, token=token), token


# send reset email test
def send_notification_email_test(person):
    return render_template('email/notification_text.html', user=person.user)


# verify token link is valid and return user id
def verify_reset_password_token(token):
    try:
        id = jwt.decode(token, app.config['SECRET_KEY'], algorithms=['HS256'])['reset_password']
    except:
        return
    return db.User.get(id=id)
