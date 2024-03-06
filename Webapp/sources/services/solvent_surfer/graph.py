import json
from typing import List

import numpy as np
import pandas as pd
from bokeh.embed import json_item
from bokeh.models import (
    Button,
    ColumnDataSource,
    CustomJS,
    Div,
    PointDrawTool,
    Range1d,
    Select,
    TapTool,
)
from bokeh.models.layouts import Tabs
from bokeh.plotting import figure

from .graph_utils import (
    CHEM21_colours,
    add_hovertool,
    colour_bar,
    control_point_colours,
    descriptor_colours,
    exit_button_setup,
    get_exp_data_dict,
    get_exp_data_fields,
    get_reaction_classes,
    interactive_mode_button_setup,
    reset_button_setup,
    save_button_setup,
    scatter_plot,
    tab_setup,
    yield_colours,
)
from .javascript_callbacks import (
    get_class_callback,
    get_click_point_callback,
    get_colour_callback,
    get_name_callback,
    get_update_callback,
)


def editable_bokeh_graph(
    df: pd.DataFrame,
    names: pd.Series,
    descriptors: List,
    colour_name: str,
    mode: str,
    control_points: np.array,
    point: List,
    r_class="",
    name="",
    return_json=True,
):
    """
    this function plots the interactive graph for the solvent surfer using bokeh. It relies on functions
    in .graph_utils.py and javascript_callbacks.py

    Args:
        df: pd.DataFrame, dataframe of the kPCA data to plot
        names: pd.Series, contains all names of solvents in the solvent surfer
        descriptors: List, contains all descriptors from initial dataset, not just reaction specific descriptors
        colour_name: str, descriptor to use to colour datapoints by
        mode: str, mode to plot the graph in ('start_up' or 'interactive' to affect graph plotting)
        control_points: np.array, array of 0's with 1's where control points are (this is added to df as a new column)
        point: List, index of point selected from solvent surfer, for updating plots on click point event
        r_class="", reaction class that appears in class_dropdown. Defaults to empty string
        name="", name of solvent that appears in name_dropdown. Defaults to empty string
        return_json=True, Bool, determines whether function returns JSON object or not. Defaults to True
    """
    # organise dataframe
    if "interactive" in mode:
        df["control_point"] = control_points

    df["names"] = names

    df = df.replace("", "NA")
    df = ColumnDataSource(
        data=df.drop([x for x in df.columns if x == "level_0"], axis=1)
    )

    # set up plot with dimensions
    p1 = figure(tools="tap", sizing_mode="stretch_width", height=500)

    # add the hover tool to show names
    add_hovertool(p1)

    # set up values for dropdowns
    names_list = names.tolist()

    names_list.insert(0, "")

    descriptors.insert(0, "CHEM21") if "CHEM21" not in descriptors else descriptors

    exp_data = get_exp_data_dict()

    exp_data_fields = get_exp_data_fields()

    if r_class in exp_data:
        for exp in exp_data[r_class]:
            descriptors.insert(1, exp) if exp not in descriptors else descriptors

    # set up dropdowns using values above
    dropdown_class = Select(
        options=get_reaction_classes(),
        value=r_class,
        sizing_mode="stretch_width",
    )

    dropdown_name = Select(options=names_list, value=name, sizing_mode="stretch_width")

    dropdown_colour = Select(
        options=descriptors, value=colour_name, sizing_mode="stretch_width"
    )

    # label dropdowns
    class_label = Div(text="""<br><b>Reaction Class</b><br>""")
    name_label = Div(text="""<br><b>Substitution Target</b><br>""")
    colour_label = Div(text="""<br><b>Change Colour</b><br>""")

    # colour control points red if applicable
    mapper_control_point = control_point_colours()

    # define default colour
    mapper = CHEM21_colours()

    # update colour mapper depending on colour_name
    if colour_name in exp_data_fields and colour_name in exp_data[r_class]:
        mapper = yield_colours(colour_name)
        # create and add colour bar
        color_bar = colour_bar(colour_name, mapper)
        p1.add_layout(color_bar, "right")

    elif colour_name != "CHEM21" and colour_name in descriptors:
        mapper = descriptor_colours(colour_name, df)
        # create and add colour bar
        color_bar = colour_bar(colour_name, mapper)
        p1.add_layout(color_bar, "right")

    else:
        colour_name = 'CHEM21'

    # set up scatter plot
    scatter = scatter_plot(p1, df, mapper_control_point, mapper, colour_name)
    p1.xaxis.axis_label = "PC1"
    p1.yaxis.axis_label = "PC2"

    # fixed range (for dragging points), 10 % bigger than max/min values
    p1.x_range = Range1d(
        min(df.data["PC1"]) - np.abs(min(df.data["PC1"]) * 0.1),
        max(df.data["PC1"]) + np.abs(min(df.data["PC1"]) * 0.1),
    )

    p1.y_range = Range1d(
        min(df.data["PC2"]) - np.abs(min(df.data["PC2"]) * 0.1),
        max(df.data["PC2"]) + np.abs(min(df.data["PC2"]) * 0.1),
    )

    # define buttons and tools
    interactive_mode = interactive_mode_button_setup(
        df, control_points, dropdown_class.value, dropdown_colour.value
    )

    reset = reset_button_setup(dropdown_class.value, dropdown_colour.value)

    exit_interactive = exit_button_setup(
        df, names, name, descriptors, dropdown_class.value, dropdown_colour.value
    )

    save_graph = save_button_setup(
        df, control_points, descriptors, dropdown_class.value, dropdown_colour.value
    )

    tap_tool = p1.select(type=TapTool)
    draw_tool = PointDrawTool(renderers=[scatter], empty_value="black", add=False)

    # if not in interactive mode, only view interactive_mode button
    buttons = interactive_mode

    # define callbacks on dropdown changes
    dropdown_colour.js_on_change(
        "value", get_colour_callback(df, dropdown_name.value, dropdown_class.value)
    )

    dropdown_class.js_on_change(
        "value",
        get_class_callback(df, point, dropdown_name.value, dropdown_colour.value),
    )

    dropdown_name.js_on_change(
        "value",
        get_name_callback(df, names, dropdown_class.value, dropdown_colour.value),
    )

    # define callback for clicking points
    tap_tool.callback = get_click_point_callback(
        df, dropdown_name.value, dropdown_class.value, dropdown_colour.value
    )

    # define callback for updating PCA once points are dragged
    df.js_on_change(
        "data",
        get_update_callback(
            df, names, name, descriptors, dropdown_class.value, dropdown_colour.value
        ),
    )

    df.js_on_change(
        "data",
        CustomJS(args=dict(draw=draw_tool), code="draw.active = !draw.active"),
    )

    df.selected.indices = point

    # require reaction class to be selected
    if mode == "start_up":
        start_up_tab = tab_setup(
            p1,
            mode,
            class_label,
            dropdown_class,
            name_label,
            dropdown_name,
            colour_label,
            dropdown_colour,
            buttons,
        )

        tabs = Tabs(tabs=[start_up_tab])
        return json_item(tabs)

    if "interactive" in mode:
        # disable interactive_mode button if in interactive mode
        interactive_mode = Button(label="Enter interactive mode", disabled=True)

        # activate draw tool for click and drag
        p1.add_tools(draw_tool)
        p1.toolbar.active_tap = draw_tool

        # add reset, exit and save buttons when in interactive mode
        buttons = [interactive_mode, reset, exit_interactive, save_graph]

    graph_tab = tab_setup(
        p1,
        mode,
        class_label,
        dropdown_class,
        name_label,
        dropdown_name,
        colour_label,
        dropdown_colour,
        buttons,
    )

    tabs = Tabs(tabs=[graph_tab])

    if not return_json:
        return json_item(tabs)

    else:
        return json.dumps(json_item(tabs))
