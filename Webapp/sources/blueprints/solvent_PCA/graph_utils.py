from typing import Dict, List, Union

import numpy as np
import pandas as pd
from bokeh.events import ButtonClick
from bokeh.layouts import column, row
from bokeh.models import (
    Button,
    CategoricalColorMapper,
    ColorBar,
    Div,
    HoverTool,
    Select,
    TabPanel,
)
from bokeh.palettes import Viridis256
from bokeh.plotting import figure
from bokeh.transform import linear_cmap

from .javascript_callbacks import (
    get_exit_interactive_callback,
    get_interactive_setup_callback,
    get_reset_interactive_callback,
    get_save_callback,
)


def get_html_dict() -> Dict:
    """
    Gets a python dictionary with keys of reaction class and values of .html files for the 'Reaction Class' tab in the
    solvent surfer.

    Returns:
        html_dict: Dict, containing .html files for each reaction class
    """

    html_dict = {
        "Amide bond formation": "reaction_class/amide_bond_formation.html",
        "Suzuki-Miyaura coupling": "reaction_class/suzuki.html",
        "Baylis-Hillman": "reaction_class/baylis.html",
        "Buchwald-Hartwig coupling": "reaction_class/buchwald.html",
        "Heck cross-coupling": "reaction_class/heck.html",
        "Grignard": "reaction_class/grignard.html",
        "SNAr/SN2": "reaction_class/snar.html",
        "Alkene metathesis": "reaction_class/metathesis.html",
        "Ester hydrolysis": "reaction_class/ester.html",
        "Alcohol oxidation": "reaction_class/alcohol_oxidation.html",
        "Other": "reaction_class/other.html",
    }

    return html_dict


def get_exp_data_dict() -> Dict:
    """
    Gets python dictionary with keys of reaction class and values of column names for yield data for that reaction class.

    Returns:
        exp_data: Dict, containing column names for yield data organised by reaction class
    """
    exp_data = {
        "Grignard": ["Grignard"],
        "Amide bond formation": ["Amide Coupling 2", "Amide Coupling 1"],
        "Heck cross-coupling": ["Heck"],
        "Baylis-Hillman": ["Baylis-Hillman"],
        "Suzuki-Miyaura coupling": ["Suzuki-Miyaura Ni", "Suzuki-Miyaura Pd"],
        "Alkene metathesis": ["Alkene Metathesis"],
        "SNAr/SN2": ["SNAr"],
        "Buchwald-Hartwig coupling": ["Buchwald-Hartwig"],
    }

    return exp_data


def get_exp_data_fields() -> List:
    """
    Gets a list of column names for all yield data included in the solvent surfer data

    Returns:
        exp_data_fields: List, column names for all yield data in solvent surfer
    """
    exp_data_fields = [
        "Grignard",
        "Amide Coupling 1",
        "Amide Coupling 2",
        "Heck",
        "Baylis-Hillman",
        "Suzuki-Miyaura Pd",
        "Suzuki-Miyaura Ni",
        "Buchwald-Hartwig",
        "Alkene Metathesis",
        "SNAr",
    ]

    return exp_data_fields


def get_reaction_classes() -> List:
    """
    Gets a list of reaction classes included in solvent surfer

    Returns:
        List, contains all reaction classes contained in solvent surfer
    """
    reaction_classes = [
        "",
        "Amide bond formation",
        "Alcohol oxidation",
        "Alkene metathesis",
        "Baylis-Hillman",
        "Buchwald-Hartwig coupling",
        "Ester hydrolysis",
        "Grignard",
        "Heck cross-coupling",
        "SNAr/SN2",
        "Suzuki-Miyaura coupling",
        "Other",
    ]

    return reaction_classes


def CHEM21_colours() -> Dict:
    """
    This function gets the bokeh colour mapper for the CHEM21 colour.

    Returns:
        Dict, contains field name and colour map for plotting datapoints
    """
    mapper = CategoricalColorMapper(
        factors=[
            "No Ranking",
            "Recommended",
            "Problematic",
            "Hazardous",
            "Highly Hazardous",
        ],
        palette=["#808080", "#008000", "#ffd700", "#FF0000", "#8b0000"],
    )

    return {"field": "CHEM21", "transform": mapper}


def yield_colours(colour_name: str) -> linear_cmap:
    """
    Gets colour map to colour solvent surfer points according to yield data from range 0-100.

    Args:
       colour_name: str, name of yield data to colour points by (can be selected from colour_dropdown)

    Returns:
        mapper: linear_cmap, mapper Object for colouring solvent surfer datapoints
    """
    mapper = linear_cmap(field_name=colour_name, palette=Viridis256, low=0, high=100)

    return mapper


def descriptor_colours(colour_name: str, df: pd.DataFrame) -> linear_cmap:
    """
    gets colour map for solvent surfer graph to colour data points by descriptors that are not CHEM21 or yield data.

    Args:
        colour_name: str, any descriptor included in solvent surfer data, but NOT 'CHEM21' or yield data columns
        df: pd.DataFrame, contains PCA data and all descriptors

    Returns:
        mapper: linear_cmap, mapper Object for colouring solvent surfer datapoints
    """
    mapper = linear_cmap(
        field_name=colour_name,
        palette=Viridis256,
        low=min([x for x in df.data[colour_name] if not isinstance(x, str)]),
        high=max([x for x in df.data[colour_name] if not isinstance(x, str)]),
    )

    return mapper


def control_point_colours() -> linear_cmap:
    """
    Gets colour map that colours points red if they are a control point (value of 1 in control_point column)

    Returns:
        cp_colour: linear_cmap, colour map for control points
    """
    cp_colour = linear_cmap(
        field_name="control_point", palette=["rgba(0, 0, 0, 0)", "red"], low=0, high=1
    )

    return cp_colour


def colour_bar(colour_name: str, mapper: linear_cmap) -> ColorBar:
    """
    Returns a colour bar for solvent surfer plot if colour_name is not 'CHEM21'

    Args:
        colour_name: str, colour name to label colour bar with (cannot be CHEM21)
        mapper: linear_cmap, colour map for bokeh plot

    Returns:
        bar: ColorBar, colour bar that is shown on solvent surfer
    """
    bar = ColorBar(
        color_mapper=mapper["transform"],
        width=8,
        location=(0, 0),
        title=colour_name,
    )

    return bar


def add_hovertool(fig: figure):
    """
    Adds hover tool to solvent surfer to show cost, melting point and boiling point when a solvent is hovered over

    Args:
        fig: bokeh.plotting.figure Object, initial bokeh figure to add hover tool to
    """
    fig.add_tools(
        HoverTool(
            tooltips=[
                ("Solvent", "@names"),
                ("Cost", "@cost"),
                ("MP", "@MP{0.0}"),
                ("BP", "@BP{0.0}"),
            ],
            mode="mouse",
        )
    )


def scatter_plot(
    fig: figure,
    source: pd.DataFrame,
    line_colour: linear_cmap,
    mapper: linear_cmap,
    colour: str,
) -> figure.scatter:
    """
    plots PCA data onto bokeh figure using right colour mappers

    Args:
        fig: bokeh.plotting.figure Object, initial bokeh figure to add datapoints to
        source: pd.DataFrame, contains all PCA data for plotting
        line_colour: linear_cmap, output of .control_point_colours()
        mapper: linear_cmap, colour mapper for chosen colour
        colour: str, colour chosen from colour_dropdown

    Returns:
        scatter: bokeh.plotting.figure.scatter, scatter plot containing kPCa data
    """

    if colour == "CHEM21":
        scatter = fig.scatter(
            "PC1",
            "PC2",
            source=source,
            alpha=0.7,
            line_alpha=1,
            line_color=line_colour,
            legend_group="CHEM21",
            line_width=4,
            color=mapper,
            size=14,
        )
        fig.add_layout(fig.legend[0], "right")

    else:
        scatter = fig.scatter(
            "PC1",
            "PC2",
            source=source,
            alpha=0.7,
            line_alpha=1,
            line_color=line_colour,
            line_width=4,
            color=mapper,
            size=14,
        )

    fig.toolbar_location = None

    return scatter


def interactive_mode_button_setup(
    df: pd.DataFrame,
    control_points: np.array,
    dropdown_class: str,
    dropdown_colour: str,
) -> Button:
    """
    Sets up solvent surfer button to enter interactive mode by linking to javascript callback

    Args:
        df: pd.DataFrame, dataframe containing PCA data from plot
        control_points: np.array, binary array with len(df) where 1 shows control points and all other points are 0
        dropdown_class: str, value selected from the class dropdown
        dropdown_colour: str, value selected from the colour dropdown

    Returns:
        interactive_mode: Button, button that calls interactive_setup callback on click
    """
    interactive_mode = Button(label="Enter interactive mode")

    interactive_mode.js_on_event(
        ButtonClick,
        get_interactive_setup_callback(
            df, control_points, dropdown_class, dropdown_colour
        ),
    )

    return interactive_mode


def reset_button_setup(dropdown_class: str, dropdown_colour: str) -> Button:
    """
    Sets up solvent surfer button to reset interactive mode by linking to javascript callback

    Args:
        dropdown_class: str, value selected from the class dropdown
        dropdown_colour: str, value selected from the colour dropdown

    Returns:
        reset: Button, button that calls reset_interactive callback on click
    """
    reset = Button(label="Reset")

    reset.js_on_event(
        ButtonClick, get_reset_interactive_callback(dropdown_class, dropdown_colour)
    )

    return reset


def exit_button_setup(
    df: pd.DataFrame,
    names: pd.Series,
    name: str,
    descriptors: List,
    dropdown_class: str,
    dropdown_colour: str,
) -> Button:
    """
    Sets up solvent surfer button to exit interactive mode by linking to javascript callback

    Args:
        df: pd.DataFrame, dataframe containing PCA data from plot
        names: pd.Series, contains all solvent names from solvent surfer
        name: str, name selected from name dropdown
        descriptors: List, all descriptors contained in the original dataset
        dropdown_class: str, value selected from the class dropdown
        dropdown_colour: str, value selected from the colour dropdown

    Returns:
        interactive_mode: Button, button that calls exit_interactive callback on click
    """
    exit_interactive = Button(label="Exit")

    exit_interactive.js_on_event(
        ButtonClick,
        get_exit_interactive_callback(
            df, names, name, descriptors, dropdown_class, dropdown_colour
        ),
    )

    return exit_interactive


def save_button_setup(df, control_points, descriptors, dropdown_class, dropdown_colour):
    """
    Sets up solvent surfer button to save edited graphs by linking to javascript callback

    Args:
        df: pd.DataFrame, dataframe containing PCA data from plot
        control_points: np.array, binary array with len(df) where 1 shows control points and all other points are 0
        descriptors: List, all descriptors contained in the original dataset
        dropdown_class: str, value selected from the class dropdown
        dropdown_colour: str, value selected from the colour dropdown

    Returns:
        save_graph: Button, button that calls save_graph callback on click
    """
    save_graph = Button(label="Save")

    save_graph.js_on_event(
        ButtonClick,
        get_save_callback(
            df, descriptors, control_points, dropdown_class, dropdown_colour
        ),
    )

    return save_graph


def tab_setup(
    fig: figure,
    mode: str,
    class_label: Div,
    dropdown_class: Select,
    name_label: Div,
    dropdown_name: Select,
    colour_label: Div,
    dropdown_colour: Select,
    buttons: List,
) -> TabPanel:
    """
    Sets up tab panel for PCA graph

    Args:
        fig: figure, bokeh figure to show on tab
        mode: str, mode graph was plotted in
        class_label: Div, label for class dropdown
        dropdown_class: Select, class dropdown
        name_label: Div, label for name dropdown
        dropdown_name: Select, name dropdown
        colour_label: Div, label for colour dropdown
        dropdown_colour: Select, colour dropdown
        buttons: List, list of buttons to show below graph
    """
    if mode == "start_up":
        tab1 = TabPanel(child=row(column(class_label, dropdown_class)), title="PCA")

    else:
        tab1 = TabPanel(
            child=column(
                row(
                    column(class_label, dropdown_class),
                    column(name_label, dropdown_name),
                    column(colour_label, dropdown_colour),
                    sizing_mode="stretch_width",
                ),
                fig,
                row(buttons),
                sizing_mode="stretch_width",
            ),
            title="PCA",
        )

    return tab1
