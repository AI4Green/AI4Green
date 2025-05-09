import dash
from flask import render_template
from markupsafe import Markup


class Dash(dash.Dash):
    """Edit the Dash class to integrate with base HTML"""

    def interpolate_index(
        self,
        metas="",
        title="",
        css="",
        config="",
        scripts="",
        app_entry="",
        favicon="",
        renderer="",
    ):
        return render_template(
            "retrosynthesis/dash.html",
            metas=Markup(metas),
            css=Markup(css),
            dash_config=Markup(config),
            scripts=Markup(scripts),
            app_entry=Markup(app_entry),
            renderer=Markup(renderer),
        )
