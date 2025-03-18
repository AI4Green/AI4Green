import json
import os
from typing import Optional

import pandas as pd
import plotly
import plotly.express as px
from flask import Response, render_template
from flask_login import login_required
from sources.auxiliary import get_notification_number, get_workgroups

from . import solvent_guide_bp


def get_radar_plot(s: str, h: str, e: str) -> str:
    df = pd.DataFrame(dict(r=[s, h, e], theta=["S", "H", "E"]))
    fig = px.line_polar(df, r="r", theta="theta", line_close=True, range_r=[0, 10])
    fig.update_traces(fill="toself")
    fig.update_xaxes(range=[0, 5])
    fig.update_layout(
        margin=dict(l=20, r=20, t=40, b=20),
    )
    fig.update_layout(height=250, width=250)
    return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)


@solvent_guide_bp.route("/solvent_guide", methods=["GET", "POST"])
@solvent_guide_bp.route("/solvent_guide/<sol>", methods=["GET", "POST"])
@login_required
@solvent_guide_bp.doc(security="sessionAuth")
def solvent_guide(sol: Optional[str] = None) -> Response:
    """
    Renders the solvent guide page. Includes a solvent card if one is provided as an argument

    Args:
        sol: solvent name

    Returns:
        flask.Response: renders the solvent guide template
    """
    CHEM21 = pd.read_csv(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "CHEM21_full.csv")
    )
    CHEM21 = CHEM21.sort_values(by="Family")
    CHEM21 = CHEM21.sort_values(by="Solvent")
    CHEM21 = CHEM21.fillna("")
    solvents = CHEM21.to_dict("index")
    solvents = list(solvents.values())
    families = sorted(list(set(CHEM21["Family"].tolist())))
    # if from reaction table and a solvent was selected. Alternative name is name used in compound/solvent database tables
    if sol:
        if not CHEM21["Solvent Alternative Name"].str.lower().eq(sol.lower()).any():
            sol = None
        else:
            sol = CHEM21[CHEM21["Solvent Alternative Name"].str.lower() == sol.lower()][
                "Number"
            ].iloc[0]
    return render_template(
        "solvent_guide.html",
        solvents=solvents,
        families=families,
        sol=sol,
    )


@solvent_guide_bp.route("/solvent_guide_help", methods=["GET", "POST"])
@login_required
@solvent_guide_bp.doc(security="sessionAuth")
def solvent_guide_help() -> Response:
    """
    Renders the solvent guide help page

    Returns:
        flask.Response: renders the solvent guide help template with information about the solvent guide
    """
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    return render_template(
        "solvent_guide_help.html",
        workgroups=workgroups,
        notification_number=notification_number,
        graphJSON=[get_radar_plot(8, 3, 5)],
    )
