import ast
import os
from typing import Dict, List

import numpy as np
import pandas as pd
import json
from flask import Response, render_template, request
from sources import models
from sources.extensions import db

from .best_closest import get_closest
from .interactive.interactive_embeddings import GetEmbedding
from .interactive.Embedder import cPCA
from sources.blueprints.main.routes import get_custom_colours


def extract_names(name_list: List) -> List:
    """
    adds blank string to match name_dropdown in solvent surfer
    """
    names = [""]
    names.extend(name_list)
    names = [x for x in names]

    return names


def get_surfer_names() -> (List, List, List):
    """
    Gets solvent names and alternative solvent names from solvent data

    Returns:
        names: List, primary names for each solvent
        alt_names: List, alternative names for each solvent
        all_names: List, combination of names and alt_names
    """
    # get names and alternative names from the solvent surfer .csv file
    solvent_data = pd.read_csv(
        os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "data", "solvent_data_PCA.csv"
        )
    )

    solvent_data = solvent_data.fillna('')
    names = extract_names(solvent_data["Solvent"].values.tolist())
    alt_names = extract_names(solvent_data["Alternative Name"].values.tolist())

    all_names = names + alt_names

    return names, alt_names, all_names


def find_point_from_name(name: str, names: List) -> List:
    """
    finds index of name in names and sets this index as point

    Args:
        name: str, name to get point for
        names: List, list of primary names os solvent sin solvent surfer

    Returns:
        point: List, index of name in names
    """
    upper_names = [x.upper() for x in names]

    point = [upper_names.index(name.upper()) - 1]

    return point


def find_reaction_table_mode_variables(
    r_class: str, name: str, names: List, alt_names: List
) -> (str, str, List):
    """
    Organises graph variables to set up solvent surfer from reaction table

    Args:
        r_class: str, reaction class selected in dropdown in reaction table
        name: str, solvent name selected from reaction table
        names: List, all name included in solvent surfer
        alt_names: List, alternative solvent names

    Returns:
        r_class: same as input unless r_class='-select-'
        name: solvent name
        point: index of name in names to label point in surfer
    """
    if r_class == "-select-":
        r_class = "Other"

    upper_names = [x.upper() for x in names]
    upper_alt = [y.upper() for y in alt_names]
    upper_name = name.upper()

    if upper_name in upper_names:
        if name == '':
            point = []
        else:
            point = [upper_names.index(upper_name) - 1]
            name = names[point[0] + 1]

    elif upper_name in upper_alt:
        point = [upper_alt.index(upper_name) - 1]
        name = names[point[0] + 1]

    else:
        point = []

    return r_class, name, point


def sort_mode_variables(
    mode: str, name: str, names: List, alt_names: List, point: List, r_class: str
) -> (str, str, List):
    """
    Sorts graph variables depending on mode

    Args:
        mode: str, mode to get variables for
        name: str, name selected from point click or name dropdown
        names: List, all solvent names included in solvent surfer
        alt_names: List, alternative solvent names
        point: List, index of name selected in names list
        r_class: str, reaction class selected in dropdown

    Returns:
        r_class: str, updated r_class if changed, or input value if not
        name: str, name extracted from names list using point
        point: List, name index in names for selecting point
    """

    if "name" in mode:
        point = find_point_from_name(name, names)

    elif mode == "point_change":
        name = names[point[0] + 1]

    elif mode == "from_reaction_table":
        r_class, name, point = find_reaction_table_mode_variables(
            r_class, name, names, alt_names
        )

    return r_class, name, point


def get_mode(
    mode: str,
) -> (str, List, str, str, List, List, str, pd.DataFrame, np.array):
    """
    Gets data from front end and arranges variables depending on mode.
    Relies on all functions in this file that organise mode variables

    Args:
        mode: str, mode with which to organise variables

    Returns:
        colour_name: str, name of colour selected in from drop down (if any)
        point: List, point selected from front end (if any)
        r_class: str, reaction class selected from front end (if any)
        name: str, name selected from front end (if any)
        names: List, contains all names of solvents included in surfer
        mode: str, mode with which to organise variables
        df: pd.DataFrame, data used to plot data
        control_points: np.array, binary array of len(df) where 1's denote control points
    """
    # organises variables for the PCA graph depending on the mode
    names, alt_names, all_names = get_surfer_names()

    # modes for get graph
    colour_name = "CHEM21"
    point = []
    r_class = ""
    name = ""
    df = pd.DataFrame()
    control_points = {}

    mode_dict = {
        "colour_selected": colour_name,
        "class_selected": r_class,
        "name_selected": name,
        "point": point,
        "data": df,
        "control_points": control_points,
        "graph_id": "graph_id",
    }

    if mode != "start_up":
        data = request.get_json()
        for key, val in mode_dict.items():
            try:
                var = data[key]

                mode_dict[key] = var

            except KeyError:
                pass

        (
            colour_name,
            r_class,
            name,
            point,
            df,
            control_points,
            graph_id,
        ) = mode_dict.values()

        df = bokeh_to_pandas(df)

        r_class, name, point = sort_mode_variables(
            mode, name, names, alt_names, point, r_class
        )

    if r_class == "":
        mode = "start_up"
        point = []
        name = ""

    return colour_name, point, r_class, name, names, mode, df, control_points


def get_update_interactive_mode(
    mode: str,
) -> (pd.DataFrame, str, str, List, np.array, cPCA, bool):
    """
    Organises graph variables for interactive mode in solvent surfer.
    Either gets variables from the database (for on_load) or from the front end (for update interactive)

    Args:
        mode: str, on_load or update interactive

    Returns:
        new_data: pd.DataFrame, updated data from graph or database
        r_class: str, reaction class selected
        colour: str, colour_selected
        descriptors: List, all descriptors included in original dataset
        control_points: nmp.array, binary array of len(new data) where 1's denote control points
        interactive_data: cPCA object, for updating PCA if more points are moved
        return_json: Bool, if True returns PCA graph as a JSON object. defaults as True
    """
    data = request.get_json()

    control_points = {}
    interactive_data = []
    embedding_dict = {}

    if mode == "on_load":
        query = db.session.query(models.PCAGraph).filter(
            models.PCAGraph.id == data["graph_id"]
        )
        r_class = query[0].r_class
        colour = query[0].colour_selected

        cp = json.loads(query[0].control_points)

        embedding_dict = json_decoder(json.loads(query[0].embedding_algorithm))

        # convert strings to int after reading json
        control_points = {int(k): v for k, v in cp.items()}

        new_data = json.loads(query[0].graph_data)
        descriptors = ast.literal_eval(query[0].descriptors)

        return_json = False

    else:
        new_data = data["data"]
        r_class = data["r_class"]
        colour = data["colour_selected"]
        descriptors = data["descriptors"]
        return_json = True

    return (
        new_data,
        r_class,
        colour,
        descriptors,
        control_points,
        interactive_data,
        embedding_dict,
        return_json,
    )


def get_suggest_solvent_table(point: List, df: pd.DataFrame) -> (Response, List):
    """
    Renders html template of suggested_solvent_table.

    Args:
        point: List, point selected for potential substitution
        df: pd.DataFrame, data used for PCA plot

    Returns:
        suggest_solvent_table: Response, rendered html template of closest solvent to selected point
        best_solvents: List, solvents arranged by Euclidean distance to point
    """
    # renders suggest solvent table template based on the selected point. needs output df of get_pca()

    if len(point) > 0:
        colour_dict = get_custom_colours().get_json()['colours']
        best_df = get_closest(point[0], df)
        best_solvents = best_df.to_dict("index")
        best_solvents = [value for value in best_solvents.values()]
        try:
            sol = best_df.iloc[0]["Solvent"]
        except KeyError:
            sol = best_df.iloc[0]["names"]

        # this is the solvent selected
        suggest_solvent_table = render_template(
            "_suggest_solvent_table.html", solvents=best_solvents, sol=sol, colours=colour_dict
        )
    else:
        suggest_solvent_table = None
        best_solvents = []

    return suggest_solvent_table, best_solvents


def search_suggest_solvent_table(best_solvents: List) -> (List, bool):
    """
    searches best_solvents for close by solvent that are more sustainable than selected solvent.

    Args:
        best_solvents: List, closest solvents to selected point arranged by Euclidean distance

    Returns:
        alternatives: List, more sustainable solvents close by to selected solvent (if any)
        substitutions: bool, sets to True if sustainable solvent switch is found
    """
    # searches suggested solvent table for sustainable switches and returns suggested substitutions

    rankings = {
        "Recommended": 1,
        "Problematic": 2,
        "Hazardous": 3,
        "Highly Hazardous": 4,
    }
    alternatives = []
    substitutions = False
    chosen = best_solvents[0]

    if chosen["CHEM21"] != "Recommended":
        for i in best_solvents:
            if 0 < i["distances"] < 0.06:
                if i["CHEM21"] != "No Ranking":
                    if rankings[i["CHEM21"]] < rankings[chosen["CHEM21"]]:
                        alternatives.append(i["Solvent"])

    if len(alternatives) > 0:
        substitutions = True

    return alternatives, substitutions


def bokeh_to_pandas(bokeh_data: Dict) -> pd.DataFrame:
    """
    this function converts a bokeh data object to a pandas dataframe.
    """
    pandas_df = pd.DataFrame()

    for col in bokeh_data.keys():
        val = bokeh_data[col]
        if not isinstance(val, list):
            vals = [*val.values()]
            for _ in range(3):
                del vals[-1]
            pandas_df[col] = vals
        else:
            pandas_df[col] = val

    return pandas_df


def update_kPCA(
    kPCA: GetEmbedding, control_points_kPCA: Dict
) -> (GetEmbedding, np.array):
    """
    updates kPCA object on control_point change

    Args:
        kPCA: GetEmbedding Object, defined in routes
        control_points_kPCA: control point dictionary from routes

    Returns:
        kPCA: GetEmbedding Object, updated object
        points: np.array, new data points for replotting graph
    """
    kPCA.embedding_algorithm.update_control_points(control_points_kPCA)
    kPCA.embedding_algorithm.finished_relocating()
    points = kPCA.update()

    return kPCA, points


def control_points_for_df(control_points: Dict, names: List) -> np.array:
    """
    constructs a binary array of len(names) to identify control points. can be added to the graph dataframe to identify control points.
    this is used to generate the colour maps for control points and for other identification methods.

    Args:
        control_points: Dict, dictionary of control points defined in routes.
        names: List, contains all solvent names in solvent surfer
    """
    c_points = np.zeros(len(names))
    for keys in control_points.keys():
        c_points[keys] = 1

    return c_points

class NumpyEncoder(json.JSONEncoder):
    """
    Class for converting np.arrays to lists for JSON sterilization.
    Used only in jsonify_embedding_params
    """
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, set):
            return list(obj)
        return json.JSONEncoder.default(self, obj)


def json_decoder(embedding_dict: Dict) -> Dict:

    """
    Converts the embedding params from the database (which are saved as JSON)
    to the appropriate types.
    Args:
        embedding_dict: Dict, JSON dict that is saved to the database as part of PCA_graph model

    Returns:
        decoded_dict: Dict, dictionary containing information to set up embedding with appropriate types
    """

    decoded_dict = {}

    array_keys = ["data", "X", "Y", "projection_matrix", "projection", "pca_projection"]
    tuple_keys = ["kernel_sys", "quad_eig_sys", "quad_eig_sys_original"]

    for key, val in embedding_dict.items():
        if key in array_keys:
            decoded_dict[key] = np.array(val)

        elif key in tuple_keys:
            if key != "quad_eig_sys":
                bro = tuple(np.array(x) for x in val)

            else:
                bro = (np.array(val[0]), [np.array(x) for x in val[1]])

            decoded_dict[key] = bro

        else:
            decoded_dict[key] = val

    return decoded_dict


def jsonify_embedding_params(embedding_algorithm: cPCA) -> json.dumps:
    """
    Extracts params from cPCA object used to generate the embedding and
    converts to JSON to save to the database

    Args:
        embedding_algorithm: cPCA object, the object used to generate the PCA embedding

    Returns:
         JSON object, embedding parameters sterilized as JSON. This is saved to the database
    """

    embedding_algorithm_dict = embedding_algorithm.__dict__

    # remove python classes from dict (these can't be converted to JSON)
    del embedding_algorithm_dict['parent']
    del embedding_algorithm_dict['embedder']

    return json.dumps(embedding_algorithm_dict, cls=NumpyEncoder)
