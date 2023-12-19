import ast
import json
import os
from datetime import date, datetime

import numpy as np
import pandas as pd
from flask import Response, jsonify, render_template, request
from flask_login import current_user, login_required
from sources import models
from sources.services.person import from_current_user_email
from sources.services.PCA_graph import list_user_graphs, PCAgraph_from_id
from sources.auxiliary import get_notification_number, get_workgroups
from sources.extensions import db

from . import solvent_PCA_bp
from .embedding import get_PCA, get_r_class_descriptors
from .graph import editable_bokeh_graph
from .graph_utils import get_html_dict
from .interactive.data_analysis import getPointChanged
from .utils import (
    bokeh_to_pandas,
    control_points_for_df,
    get_mode,
    get_suggest_solvent_table,
    get_surfer_names,
    get_update_interactive_mode,
    search_suggest_solvent_table,
    update_kPCA,
)


@solvent_PCA_bp.route("/solvent_PCA", methods=["GET", "POST"])
@solvent_PCA_bp.route("/solvent_PCA/<mode>", methods=["GET", "POST"])
@login_required
def solvent_PCA() -> Response:
    """
    Loads initial solvent_PCA page which allows access to solvent surfer
    """
    # user must be logged in
    workgroups = get_workgroups()
    notification_number = get_notification_number()

    return render_template(
        "solvent_PCA.html",
        workgroups=workgroups,
        notification_number=notification_number,
    )


interactive_data_kPCA = []
control_points_kPCA = {}


@solvent_PCA_bp.route("/get_graph", methods=["GET", "POST"])
@solvent_PCA_bp.route("/get_graph/<mode>", methods=["GET", "POST"])
def get_graph(mode="start_up") -> Response:
    """
    Route plots the interactive graph.

    Args:
        mode: str, links graph variables to callbacks depending on user interaction. defaults as 'start_up'

    Returns:
        Response, Json object containing plotted graph (chart) and suggested_solvent_table (if applicable)
    """
    # Deal with graph variables depending on mode
    colour_name, point, r_class, name, names, mode, mode_df, control_points = get_mode(
        mode
    )

    # get the PC dataframe to plot the graph
    pca_df, descriptors, names, kPCA = get_PCA(r_class)

    if mode == "reset_interactive":
        control_points_kPCA.clear()
        interactive_data_kPCA.clear()

    elif mode == "on_load":
        control_points_kPCA.clear()
        interactive_data_kPCA.clear()

        control_points_kPCA.update(control_points)

    interactive_data_kPCA.append(pca_df)
    interactive_data_kPCA.append(kPCA)

    if len(mode_df) == 0:
        df = pca_df
    else:
        df = mode_df

    suggest_solvent_table = ""

    c_points = np.zeros(len(names))

    for keys in control_points_kPCA.keys():
        c_points[keys] = 1

    # get the suggested solvent table if a point has been selected
    if len(point) > 0:
        suggest_solvent_table, best_solvents = get_suggest_solvent_table(point, df)

    # plot the graph and send
    chart = editable_bokeh_graph(
        df,
        names,
        descriptors,
        colour_name,
        mode,
        c_points,
        point,
        r_class,
        name,
        return_json=False,
    )

    return json.dumps({"editableChart": chart, "table": suggest_solvent_table})


@solvent_PCA_bp.route("/on_point_click", methods=["GET", "POST"])
@solvent_PCA_bp.route("/on_point_click/<mode>", methods=["GET", "POST"])
def on_point_click(mode="point_change") -> Response:
    """
    renders suggested solvent table upon point click
    """
    # get graph variables depending on mode
    colour_name, point, r_class, name, names, mode, mode_df, _ = get_mode(mode)

    suggest_solvent_table = ""

    # get suggested solvent table and render template
    if len(point) > 0:
        suggest_solvent_table, best_solvents = get_suggest_solvent_table(point, mode_df)

    return jsonify({"suggestSolventTable": suggest_solvent_table})


@solvent_PCA_bp.route("/on_class_change", methods=["GET", "POST"])
def on_class_change() -> Response:
    """
    Renders 'about reaction class' .html files upon class change
    """
    data = request.get_json()
    r_class = data["class_selected"]

    control_points_kPCA.clear()
    interactive_data_kPCA.clear()

    # get reaction class html templates
    class_html = get_html_dict()

    if r_class == "":
        return jsonify({"aboutSolventClass": ""})

    else:
        # render template
        about_solvent_class = render_template(class_html[r_class])
        return jsonify({"aboutSolventClass": about_solvent_class})


@solvent_PCA_bp.route("/from_reaction_table", methods=["GET", "POST"])
@solvent_PCA_bp.route("/from_reaction_table/<mode>", methods=["GET", "POST"])
def from_reaction_table(mode="from_reaction_table") -> Response:
    """
    Sets up solvent surfer with variables from the reaction table.
    Also searches suggest solvent table if mode == 'check_solvents'
    """
    # get graph variables depending on mode
    colour_name, point, r_class, name, names, _, df, _ = get_mode(
        mode="from_reaction_table"
    )

    # get PCA dataframe
    df, descriptors, pca_names, kPCA = get_PCA(r_class)

    # get reaction class html templates
    class_html = get_html_dict()

    # get all solvent names included in solvent surfer
    names, alt_names, all_names = get_surfer_names()

    # get suggested solvent table
    suggest_solvent_table = ""

    if name != "":
        suggest_solvent_table, best_solvents = get_suggest_solvent_table(point, df)

    # plot graph
    editable_chart = editable_bokeh_graph(
        df,
        pca_names,
        descriptors,
        "CHEM21",
        "reaction_table",
        np.zeros(len(names)),
        point,
        r_class,
        name,
        return_json=False,
    )

    # render class html template
    about_solvent_class = render_template(class_html[r_class])

    # return graph with highlighted point if solvent surfer is opened from rxn table
    if mode == "from_reaction_table":
        return json.dumps(
            {
                "editableChart": editable_chart,
                "suggestSolventTable": suggest_solvent_table,
                "aboutSolventClass": about_solvent_class,
            }
        )

    # else return possible sustainable substitutions for alert in rxn table
    else:
        alternatives, substitutions = search_suggest_solvent_table(best_solvents)

        return jsonify(
            {
                "solvents": [x.upper() for x in names],
                "alternatives": alternatives,
                "substitutions": substitutions,
            }
        )


@solvent_PCA_bp.route("/update_interactive_graph/<mode>", methods=["GET", "POST"])
def update_interactive_graph(mode="interactive") -> Response:
    """
    Updates solvent surfer in interactive mode if points have been dragged to new positions
    """
    # this route updates the graph on dragging points

    # load variables depending on mode
    (
        new_data,
        r_class,
        colour,
        descriptors,
        control_points,
        interactive_data,
        return_json,
    ) = get_update_interactive_mode(mode)

    # get descriptors for reaction class
    class_descriptors, _ = get_r_class_descriptors(r_class)

    if mode == "on_load":
        # reset control points and PCA if graph is loaded
        control_points_kPCA.clear()
        interactive_data_kPCA.clear()

        # update control_points and kpca with loaded data
        control_points_kPCA.update(control_points)
        interactive_data_kPCA.extend(interactive_data)

        # render about solvent class html for solvent class tab
        class_html = get_html_dict()
        about_solvent_class = render_template(class_html[r_class])

    previous = interactive_data_kPCA[0]
    names = pd.Series(new_data["names"])

    # convert data from bokeh data to pandas dataframe
    full_data = bokeh_to_pandas(new_data)

    # extract first two PCs to detect dragged points
    PCs = full_data[["PC1", "PC2"]]
    PCs = PCs.dropna()

    # compare PCs of new data with old data to get control points
    point_changed = getPointChanged(PCs, previous[["PC1", "PC2"]])

    # if control_points have changed, update kPCA object and control_points
    if point_changed:
        control_points_kPCA[point_changed[0]] = point_changed[1]

        kPCA, points = update_kPCA(interactive_data_kPCA[1], control_points_kPCA)

        interactive_data_kPCA[1] = kPCA

    else:
        points = interactive_data_kPCA[1].update()

    # get list with 1 if control point, 0 if not
    c_points = control_points_for_df(control_points_kPCA, names)

    # create dataframe with updated kpca data points
    df = pd.DataFrame(data=points.T, columns=["PC1", "PC2"])

    # merge with initial descriptor data (for the colour)
    df = pd.concat([df, full_data.drop(["PC1", "PC2"], axis=1)], axis=1)

    # replace data with new data
    interactive_data_kPCA[0] = df

    editable_chart = editable_bokeh_graph(
        df,
        names,
        descriptors,
        colour,
        "interactive",
        c_points,
        [],
        r_class,
        return_json=return_json,
    )

    # return chart and suggested_solvent_table depending on mode
    if mode != "on_load":
        return editable_chart

    else:
        return json.dumps(
            {"editableChart": editable_chart, "aboutSolventClass": about_solvent_class}
        )


@solvent_PCA_bp.route("/save_graph", methods=["GET", "POST"])
def save_graph() -> Response:
    """
    Saves edited surfers as PCA_graph model
    """
    graph_data = request.get_json()

    author = from_current_user_email()

    saved_graph = models.PCAGraph(
        graph_name=graph_data["graph_name"],
        r_class=graph_data["class_selected"],
        colour_selected=graph_data["colour_selected"],
        control_points_pkl=control_points_kPCA,
        graph_data=graph_data["graph_data"],
        descriptors=str(graph_data["descriptors"]),
        kpca_object=interactive_data_kPCA,
        time_of_creation=datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f")[:-7],
        date_of_creation=date.today(),
        creator_person=author,
        status="active",
    )

    db.session.add(saved_graph)
    db.session.commit()

    all_graphs = get_saved_graphs()

    return all_graphs


@solvent_PCA_bp.route("/get_saved_graphs", methods=["GET", "POST"])
def get_saved_graphs() -> Response:
    """
    Extracts saved graphs from the database by user id for 'Saved Graphs' tab.
    """
    user = from_current_user_email()

    graphs = list_user_graphs(user)

    graph_list = []

    for idx, graph in enumerate(graphs):
        graph_details = {
            "html_id": idx + 1,
            "graph_id": graph.id,
            "name": graph.graph_name,
            "r_class": graph.r_class,
            "colour": graph.colour_selected,
            "control_points": graph.control_points_pkl,
            "time_of_creation": graph.time_of_creation,
        }

        graph_list.append(graph_details)

    graph_details = render_template("_saved_pca_graphs.html", graph_list=graph_list)

    return (
        jsonify({"graphDetails": graph_details})
        if len(graph_list) > 0
        else jsonify("no_graphs")
    )


@solvent_PCA_bp.route("/delete_graph/<graph_id>", methods=["GET", "POST"])
def delete_graph(graph_id) -> Response:
    """
    Sets a graph's status to inactive, so it doesn't show up in saved graphs.
    Graphs are not deleted from the database in case they need to be restored.
    """
    user = from_current_user_email()

    query = PCAgraph_from_id(graph_id, user)

    query.status = "inactive"
    db.session.commit()

    return get_saved_graphs()
