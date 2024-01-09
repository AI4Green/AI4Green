from typing import List

import numpy as np
import pandas as pd
from bokeh.models import CustomJS


def class_change_callback() -> str:
    """
    String for class change callback
    """
    javascript_string = """
        var class_selected = this.value;
        fetch('/get_graph/class_change', {
        // Declare what type of data we're sending
        headers: {
            'Content-Type': 'application/json'
        },
        // Specify the method
        method: 'POST',
        // A JSON payload
        body: JSON.stringify({
        "class_selected": class_selected,
        "point": point,
        "name_selected": name,
        "colour_selected": colour_selected
            })
        })
        // we want to replot the graphs unless there are fewer than 3 descriptors (alert user in this case)
        .then(function(response) { return response.json(); })
        .then(function(item) {
        $('#pre-class').hide()
        $('#suggested_table').html(item.table).show();
        $("#editable_chart").empty();
        $('#right-panel').show();
        return Bokeh.embed.embed_item(item.editableChart, "editable_chart")
        })
        .then(function(){
            if (class_selected === ""){
                $('#pre-target').show();
                $('#pre-class').show();
            }
            else{
                $('#pre-class').hide();
            }
        })
        .then(function(){
            fetch('/on_class_change', {
                // Declare what type of data we're sending
                headers: {
                  'Content-Type': 'application/json'
                },
                // Specify the method
                method: 'POST',
                // A JSON payload
                body: JSON.stringify({"class_selected": class_selected})
            })
            .then(function(response) {
                return response.json();
            })
            .then(function(item) {
                $('#about_class').html(item.aboutSolventClass).show();
            })
        })
        """

    return javascript_string


def colour_callback() -> str:
    """
    String for colour change callback
    """
    javascript_string = """
        const data = source.data;
        const indices = source.selected.indices;
        var colour_selected = this.value;
        fetch('/get_graph/colour_change', {
        // Declare what type of data we're sending
        headers: {
          'Content-Type': 'application/json'
        },
        // Specify the method
        method: 'POST',
        // A JSON payload
        body: JSON.stringify({"colour_selected": colour_selected, "class_selected": class_selected,
                            "point": indices, "name_selected": name_selected, "data": data})
        })
        // we want to replot the graph
        .then(function(response) { return response.json(); })
        .then(function(item) {
        $("#editable_chart").empty();
        return Bokeh.embed.embed_item(item.editableChart, "editable_chart");
        })
        """

    return javascript_string


def click_point_callback() -> str:
    """
    string for click point callback
    """
    javascript_string = """
        const data = source.data;
        const indices = source.selected.indices;
        fetch('/on_point_click', {
            // Declare what type of data we're sending
            headers: {
              'Content-Type': 'application/json'
            },
            // Specify the method
            method: 'POST',
            // A JSON payload
            body: JSON.stringify({"point": indices, "class_selected": class_selected,
                                    "colour_selected": colour_selected, "data": data})
        })
        .then(function(response) {
            return response.json();
        })
        .then(function(item) {
            $('#pre-target').hide();
            $('#suggested_table').html(item.suggestSolventTable).show();
        })
        .then(function() {
            fetch('/get_graph/point_change', {
            // Declare what type of data we're sending
            headers: {
                'Content-Type': 'application/json'
            },
            // Specify the method
            method: 'POST',
            // A JSON payload
            body: JSON.stringify({"class_selected": class_selected, "colour_selected": colour_selected, "point": indices,
                                    "name_selected": name_selected, "data": data})
            })
            // we want to replot the graph
            .then(function(response) { return response.json(); })
            .then(function(item) {
            $("#editable_chart").empty();
            return Bokeh.embed.embed_item(item.editableChart, "editable_chart");
            })
        })
        """

    return javascript_string


def name_callback() -> str:
    """
    string for name change callback
    """
    javascript_string = """
        const data = source.data;
        const indices = source.selected.indices;
        var name_selected = this.value;
        fetch('/get_graph/name_change', {
        // Declare what type of data we're sending
        headers: {
            'Content-Type': 'application/json'
        },
        // Specify the method
        method: 'POST',
        // A JSON payload
        body: JSON.stringify({"name_selected": name_selected, "class_selected": class_selected,
                                "colour_selected": colour_selected, "point": indices, "data": data})
        })
        // we want to replot the graph
        .then(function(response) { return response.json(); })
        .then(function(item) {
        $("#editable_chart").empty();
        return Bokeh.embed.embed_item(item.editableChart, "editable_chart");
        })
        .then(function(){
            fetch('/on_point_click/change_name', {
                // Declare what type of data we're sending
                headers: {
                  'Content-Type': 'application/json'
                },
                // Specify the method
                method: 'POST',
                // A JSON payload
                body: JSON.stringify({"names": names, "name_selected": name_selected, "point": indices,
                                        "class_selected": class_selected, "colour_selected": colour_selected, "data": data})
            })
            .then(function(response) {
                return response.json();
            })
            .then(function(item) {
                $('#pre-target').hide();
                $('#suggested_table').html(item.suggestSolventTable).show();
            })
            .then(function(){
                if (name_selected === ""){
                    $('#suggested_table').hide();
                    $('#pre-target').show();
                }
            })
        })

        """

    return javascript_string


def interactive_setup_callback() -> str:
    """
    string for interactive setup callback
    """
    javascript_string = """
        const indices = source.selected.indices;
        const data = source.data;
        fetch('/get_graph/interactive', {
        // Declare what type of data we're sending
        headers: {
            'Content-Type': 'application/json'
        },
        // Specify the method
        method: 'POST',
        // A JSON payload
        body: JSON.stringify({"class_selected": r_class, "colour_selected": colour_selected,
                                "data": data, "control_points": control_points})
        })
        // we want to replot the graph
        .then(function(response) { return response.json(); })
        .then(function(item) {
        $('#suggested_table').hide();
        $('#pre-target').show();
        $("#editable_chart").empty();

        return Bokeh.embed.embed_item(item.editableChart, "editable_chart");
        $("#editable_chart").focus();
        })

        """
    return javascript_string


def interactive_update_callback() -> str:
    """
    String for interactive setup callback
    """
    javascript_string = """
            $(editable_chart).find("*").prop("disabled", true);
            $(editable_chart).fadeTo("fast", 0.15);
            const data = source.data;
            fetch('/update_interactive_graph/interactive', {
            // Declare what type of data we're sending
            headers: {
              'Content-Type': 'application/json'
            },
            // Specify the method
            method: 'POST',
            // A JSON payload
            body: JSON.stringify({"data": data, "r_class": r_class,
                                    "descriptors": descriptors, "names": names,
                                    "colour_selected": colour_selected})
            })
            // we want to replot the graphs unless there are fewer than 3 descriptors (alert user in this case)
            .then(function(response) { return response.json(); })
            .then(function(item) { if(typeof item === 'string') {alert(item)} else{
            $(editable_chart).empty();
            $(editable_chart).find("*").prop("disabled", false);
            $(editable_chart).fadeTo("fast", 1);
            return Bokeh.embed.embed_item(item, editable_chart)};
            })
            """
    return javascript_string


def exit_interactive_callback() -> str:
    """
    String for exit interactive callback
    """
    javascript_string = """
            const data = source.data;
            fetch('get_graph/exit', {
            // Declare what type of data we're sending
            headers: {
              'Content-Type': 'application/json'
            },
            // Specify the method
            method: 'POST',
            // A JSON payload
            body: JSON.stringify({"data": data, "class_selected": r_class, "colour_selected": colour_selected})
            })
            // we want to replot the graphs unless there are fewer than 3 descriptors (alert user in this case)
            .then(function(response) { return response.json(); })
            .then(function(item) {
            $(editable_chart).empty();
            return Bokeh.embed.embed_item(item.editableChart, editable_chart);
            })
            """
    return javascript_string


def reset_interactive_callback() -> str:
    """
    String for reset interactive callback
    """
    javascript_string = """
            fetch('get_graph/reset_interactive', {
            // Declare what type of data we're sending
            headers: {
              'Content-Type': 'application/json'
            },
            // Specify the method
            method: 'POST',
            // A JSON payload
            body: JSON.stringify({"class_selected": class_selected, "colour_selected": colour_selected})
            })
            // we want to replot the graphs unless there are fewer than 3 descriptors (alert user in this case)
            .then(function(response) { return response.json(); })
            .then(function(item) {
            $(editable_chart).empty();
            return Bokeh.embed.embed_item(item.editableChart, editable_chart);
            })
            """
    return javascript_string


def save_interactive_callback() -> str:
    """
    String for save interactive callback
    """
    javascript_string = """

    let graphName = prompt("Enter name for graph:", "Enter graph name here");
    let graphData = graph_data.data

    fetch('/save_graph', {
    // Declare what type of data we're sending
    headers: {
      'Content-Type': 'application/json'
    },
    // Specify the method
    method: 'POST',

    // A JSON payload
    body: JSON.stringify({
        "class_selected": class_selected,
        "colour_selected": colour_selected,
        "control_points": control_points,
        "descriptors": descriptors,
        "graph_name": graphName,
        "graph_data": graphData,
        })
    })

    .then(function(response) { return response.json(); })
    .then(function(item) {
        let $reactionSaveIndicator = $("#reaction-saved-indicator");
        $reactionSaveIndicator.text("Graph saved!");
        $reactionSaveIndicator.removeClass().addClass("reaction-save-success").fadeIn("fast");
        $("#no-graphs").hide();
        $("#graph-details").html(item.graphDetails).fadeIn();

        setTimeout($("#reaction-saved-indicator").fadeOut("slow"), 1000);
    })

    """

    return javascript_string


def get_class_callback(
    df: pd.DataFrame, point: List, name_dropdown: str, colour_dropdown: str
) -> CustomJS:
    """
    Sets up javascript callback for class change

    Args:
        df: pd.DataFrame, dataframe used to plot graph
        point: List, index of point selected or empty list
        name_dropdown: str, value from name dropdown
        colour_dropdown: str, value from colour dropdown

    Returns:
        callback_class: CustomJS, callback for class change
    """
    callback_class = CustomJS(
        args={
            "source": df,
            "point": point,
            "name": name_dropdown,
            "colour_selected": colour_dropdown,
        },
        code=class_change_callback(),
    )

    return callback_class


def get_colour_callback(
    df: pd.DataFrame, name_dropdown: str, class_dropdown: str
) -> CustomJS:
    """
    Sets up javascript callback for colour change

    Args:
        df: pd.DataFrame, dataframe used to plot graph
        name_dropdown: str, value from name dropdown
        class_dropdown: str, value from class dropdown

    Returns:
        callback_colour: CustomJS, callback for colour change
    """
    callback_colour = CustomJS(
        args={
            "source": df,
            "class_selected": class_dropdown,
            "name_selected": name_dropdown,
        },
        code=colour_callback(),
    )

    return callback_colour


def get_click_point_callback(
    df: pd.DataFrame, name_dropdown: str, class_dropdown: str, colour_dropdown: str
) -> CustomJS:
    """
    Sets up javascript callback for point click

    Args:
        df: pd.DataFrame, dataframe used to plot graph
        name_dropdown: str, value from name dropdown
        class_dropdown: str, value from class dropdown
        colour_dropdown: str, value from colour dropdown

    Returns:
        callback_click_point: CustomJS, callback for point click
    """
    callback_click_point = CustomJS(
        args={
            "source": df,
            "class_selected": class_dropdown,
            "colour_selected": colour_dropdown,
            "name_selected": name_dropdown,
        },
        code=click_point_callback(),
    )

    return callback_click_point


def get_name_callback(
    df: pd.DataFrame, names: pd.Series, class_dropdown: str, colour_dropdown: str
) -> CustomJS:
    """
    Sets up javascript callback for name change

    Args:
        df: pd.DataFrame, dataframe used to plot graph
        names: pd.Series, contains all solvent names in solvent surfer
        class_dropdown: str, value from class dropdown
        colour_dropdown: str, value from colour dropdown

    Returns:
        callback_name: CustomJS, callback for name change
    """
    callback_name = CustomJS(
        args={
            "source": df,
            "class_selected": class_dropdown,
            "colour_selected": colour_dropdown,
            "names": names.to_list(),
        },
        code=name_callback(),
    )

    return callback_name


def get_update_callback(
    df: pd.DataFrame,
    names: pd.Series,
    name: str,
    descriptors: List,
    class_dropdown: str,
    colour_dropdown: str,
) -> CustomJS:
    """
    Sets up javascript callback to update solvent surfer upon changing points

    Args:
        df: pd.DataFrame, dataframe used to plot graph
        names: pd.Series, contains all solvent names in solvent surfer
        name: str, name selected from name dropdown
        descriptors: List, list of all descriptors included in the original dataset
        class_dropdown: str, values from class dropdown
        colour_dropdown: str, value from colour dropdown

    Returns:
        callback_interactive_update: CustomJS, callback for updating interactive mode
    """
    callback_interactive_update = CustomJS(
        args={
            "source": df,
            "r_class": class_dropdown,
            "descriptors": descriptors,
            "names": names,
            "name_selected": name,
            "colour_selected": colour_dropdown,
        },
        code=interactive_update_callback(),
    )

    return callback_interactive_update


def get_interactive_setup_callback(
    df: pd.DataFrame,
    control_points: np.array,
    class_dropdown: str,
    colour_dropdown: str,
) -> CustomJS:
    """
    Sets up javascript callback to set up interactive mode

    Args:
        df: pd.DataFrame, dataframe used to plot graph
        control_points: np.array, binary list of len(df) where 1 denotes a control point
        class_dropdown: str, value from class dropdown
        colour_dropdown: str value from colour dropdown

    Returns:
        callback_class: CustomJS, callback for class change
    """
    callback_interactive_setup = CustomJS(
        args={
            "source": df,
            "r_class": class_dropdown,
            "colour_selected": colour_dropdown,
            "control_points": control_points,
        },
        code=interactive_setup_callback(),
    )

    return callback_interactive_setup


def get_exit_interactive_callback(
    df: pd.DataFrame,
    names: pd.Series,
    name: str,
    descriptors: List,
    class_dropdown: str,
    colour_dropdown: str,
):
    """
    Sets up javascript callback to set up interactive mode

    Args:
        df: pd.DataFrame, dataframe used to plot graph
        names: pd.Series, contains all solvent names in solvent surfer
        name: name selected from dropdown or point click
        class_dropdown: str, value from class dropdown
        colour_dropdown: str value from colour dropdown

    Returns:
        callback_exit_interactive: CustomJS, callback for exiting interactive mode
    """
    callback_exit_interactive = CustomJS(
        args={
            "source": df,
            "r_class": class_dropdown,
            "descriptors": descriptors,
            "names": names,
            "name_selected": name,
            "colour_selected": colour_dropdown,
        },
        code=exit_interactive_callback(),
    )

    return callback_exit_interactive


def get_reset_interactive_callback(class_dropdown: str, colour_dropdown: str):
    """
    Sets up javascript callback to reset edited solvent surfers

    Args:
        class_dropdown: str, value from class dropdown
        colour_dropdown: str value from colour dropdown

    Returns:
        callback_reset_interactive: CustomJS, callback for resetting edited solvent surfers
    """
    callback_reset_interactive = CustomJS(
        args={
            "class_selected": class_dropdown,
            "colour_selected": colour_dropdown,
        },
        code=reset_interactive_callback(),
    )

    return callback_reset_interactive


def get_save_callback(
    df: pd.DataFrame,
    descriptors: List,
    control_points: np.array,
    class_dropdown: str,
    colour_dropdown: str,
):
    """
    Sets up javascript callback to save edited solvent surfers

    Args:
        df: pd.DataFrame, dataframe used to plot graph
        descriptors: List, contains all descriptors from the original dataset
        control_points: np.array, binary list of len(df) where 1 denotes a control point
        class_dropdown: str, value from class dropdown
        colour_dropdown: str value from colour dropdown

    Returns:
        callback_class: CustomJS, callback to save edited solvent surfers
    """
    callback_save = CustomJS(
        args={
            "class_selected": class_dropdown,
            "colour_selected": colour_dropdown,
            "control_points": control_points,
            "graph_data": df,
            "descriptors": descriptors,
        },
        code=save_interactive_callback(),
    )

    return callback_save
