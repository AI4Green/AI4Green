from flask import render_template

from . import reaction_set_bp


@reaction_set_bp.route('/reaction_set')
def reaction_set():
    return render_template('reaction_set.html')