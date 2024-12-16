from flask import render_template

from . import reaction_set_bp


@reaction_set_bp.route('/reaction_set')
def reaction_set():
    return render_template('reaction_set.html')


@reaction_set_bp.route('/click_and_drag')
def click_and_drag():
    return render_template('click-and-drag.html')