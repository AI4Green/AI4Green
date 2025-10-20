import os
import re
from typing import Dict, List, Literal, Optional, Tuple, Union
from urllib.parse import quote

import dash
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
from dash import ALL, Input, Output, State, ctx, dcc, html
from flask import Flask, current_app
from rdkit import Chem
from sources import services

from . import classes, plotly_elements, style_sheets
from .predictive_chemistry import (
    conditions_api,
    cytoscape,
    dropdowns,
    retrosynthesis_api,
    saved_retrosyntheses,
    sustainability,
    tables,
    user_uploaded_routes,
    utils,
)

cyto.load_extra_layouts()


def init_dashboard(server: Flask) -> classes.Dash:
    """Called on app initialization. Creates a Plotly Dash dashboard."""
    # make the dash app instance
    dash_app = classes.Dash(
        server=server,
        routes_pathname_prefix="/retrosynthesis/",
        external_stylesheets=[
            "/static/css/pagestyle.css",
        ],
    )
    # read api keys and urls from the yaml file
    retrosynthesis_api_key = server.config["RETROSYNTHESIS_API_KEY"]
    retrosynthesis_base_url = server.config["RETROSYNTHESIS_API_URL"]

    # Make the page layout. Elements imported plotly_elements.py
    dash_app.layout = html.Div(
        id="viewport-container",
        className="ml-1",
        children=[
            dbc.Row(
                className="g-0",
                children=[
                    # narrow column for sidebar and larger column for remainder of page
                    dbc.Col(
                        plotly_elements.tabs,
                        id="retrosynthesis-sidebar-root",
                        className="col-4 collapse show",
                    ),
                    dbc.Col(
                        className="col-8",
                        children=[
                            plotly_elements.header_and_inputs,
                            plotly_elements.retro_tree,
                        ],
                    ),
                ],
            ),
            plotly_elements.save_modal,
            plotly_elements.new_reaction_modal,
            plotly_elements.data_storage_elements,
        ],
        style=style_sheets.WebElementsStyles.CONTENT_STYLE,
    )

    """
    Callbacks
    """

    @dash_app.callback(Output("smiles-input", "value"), Input("url", "pathname"))
    def load_imported_smiles(url_dash: str) -> str:
        """
        Read smiles from the url and check if a valid molecule
        Args:
            # inputs
            url_dash - the url potentially containing the smiles imported from the sketcher
        Returns:
             smiles-input - smiles string or an empty string if not valid
        """
        imported_smiles = url_dash.split("/retrosynthesis/")[-1]
        # only update if it is a valid smiles
        imported_smiles = utils.encodings_to_smiles_symbols(imported_smiles)
        m = Chem.MolFromSmiles(imported_smiles, sanitize=False)
        if m is None:
            return dash.no_update
        if imported_smiles:
            return imported_smiles
        else:
            return ""

    @dash_app.callback(
        Output("save-modal-workbook-dropdown", "options"),
        Output("save-modal-workbook-dropdown", "value"),  # dropdown
        Output("save-modal-workbook-dropdown", "disabled"),
        Output("new-reaction-workbook-dropdown", "options"),
        Output("new-reaction-workbook-dropdown", "value"),  # dropdown
        Output("new-reaction-workbook-dropdown", "disabled"),
        Output("reload-workbook-dropdown", "options"),
        Output("reload-workbook-dropdown", "value"),  # dropdown
        Output("reload-workbook-dropdown", "disabled"),
        Output("save-functionality-status", "data"),
        Input("url", "pathname"),
    )
    def load_user_workbooks(
        url_dash: str,
    ) -> Tuple[List[dict], int, str, List[dict], int, str, List[dict], int, str]:
        """
        Loads all workbooks the user belongs to, populating the dropdowns
        and enabling/disabling the dropdowns and saving if there are no workbooks

        Args:
            # Inputs
            url_dash: Used as input to call function upon page load

        Returns:
            A list of:
                Repeated 3 times - once for each dropdown:
                    a list of dictionaries with workbook label and id
                    initially active workbook id
                    boolean to disable or enable the dropdown
                one string to enable or disable saving
        """
        (
            workbooks_dropdown,
            initial_dropdown_value,
        ) = dropdowns.make_workbooks_dropdown_options()
        if workbooks_dropdown:
            return 3 * (workbooks_dropdown, initial_dropdown_value, False) + (
                "enabled",
            )
        # if no workbooks dropdown
        return 3 * (workbooks_dropdown, None, True) + ("disabled",)

    @dash_app.callback(
        Output("open-save-modal", "disabled"),
        Output("open-new-reaction-modal", "disabled"),
        Input("save-functionality-status", "data"),
    )
    def save_features_handler(
        save_functionality_status: Union[Literal["disabled"], Literal["enabled"]],
    ) -> [bool, bool]:
        """
        Called on page load after function: 'load_user_workbooks'
        Disables or enables the save buttons

        Args:
            # Inputs
            save_functionality_status - disabled if user is not in a workbook

        Returns:
            2 bools to disable the save buttons if the user is not in a workbook or no update
        """
        if save_functionality_status == "disabled":
            return True, True
        return dash.no_update

    """Smiles input validation"""

    @dash_app.callback(
        Output("smiles-input", "pattern"),
        Input("smiles-input", "value"),
        prevent_initial_call=True,
    )
    def smiles_check(smiles_input: str) -> str:
        """
        Validate SMILES input. If invalid, modify the pattern to highlight the input field.

        Args:
            smiles_input: SMILES string for retrosynthesis.

        Returns:
            Pattern string with '_invalid' if the input is invalid.
        """
        smiles_input = (smiles_input or "").strip()
        m = Chem.MolFromSmiles(smiles_input, sanitize=False)
        if m is None or not smiles_input:
            return f"_invalid_{re.escape(smiles_input)}_invalid"
        return re.escape(smiles_input)

    @dash_app.callback(
        Output("validated-smiles", "data"),
        State("smiles-input", "value"),
        Input("smiles-input", "pattern"),
        prevent_initial_call=True,
    )
    def valid_smiles(smiles: str, smiles_regex: str) -> str:
        """
        Pass valid SMILES to the validated smiles field.

        Args:
            smiles: Target SMILES string.
            smiles_regex: String pattern to check if SMILES is valid.

        Returns:
            SMILES string if valid, otherwise no update.
        """
        if smiles and smiles_regex and "_invalid" not in smiles_regex:
            return smiles
        return dash.no_update

    """
    Retrosynthesis process
    """

    # Global store for task results
    task_results = {}

    def retrosynthesis_process_wrapper(request_url: str) -> str:
        """
        Args:
            request_url - the url with the target smiles and api key

        Returns:
            str: The ID of the job that has been started
        """
        task_id = retrosynthesis_api.retrosynthesis_api_call(
            request_url, retrosynthesis_base_url
        )
        return task_id

    @dash_app.callback(
        Output("user-message", "children", allow_duplicate=True),
        Output("computed-retrosynthesis-uuid", "data"),
        Output("interval-container", "children", allow_duplicate=True),
        Output("loading-display", "display", allow_duplicate=True),
        State("validated-smiles", "data"),
        State("smiles-input", "pattern"),
        State("enhancement-dropdown", "value"),
        State("iterations-input", "value"),
        State("time-limit-input", "value"),
        State("max-depth-input", "value"),
        Input("btn-retrosynthesis", "n_clicks"),
        prevent_initial_call=True,
    )
    def start_new_retrosynthesis(
        validated_smiles: str,
        smiles_regex: str,
        enhancement: str,
        iterations: int,
        time_limit: int,
        max_depth: int,
        n_clicks: int,
    ) -> Tuple[str, str, Optional[dcc.Interval], str]:
        """
        Called when the user clicks the retrosynthesis button.
        Starts the retrosynthesis process in a background thread.

        Args:
            n_clicks - increments when the user clicks the retrosynthesis button and calls the function
            validated_smiles - the validated SMILES string
            smiles_regex - contains 'invalid' if SMILES are not valid and prompts user to enter valid SMILES
            enhancement - the selected enhancement type from the dropdown
            iterations - the number of iterations for the search
            time_limit - the time limit for the search (in seconds)
            max_depth - the maximum depth of the search tree

        Returns:
            a message to give the user feedback
            The generated UUID for the retrosynthesis
            children for the dcc.Interval container element
        """
        if utils.smiles_not_valid(smiles_regex):
            return "Please enter a valid SMILES", "", None, "hide"

        validated_smiles = utils.encodings_to_smiles_symbols(validated_smiles)

        iterations = iterations if iterations is not None else 100
        time_limit = time_limit if time_limit is not None else 60
        max_depth = max_depth if max_depth is not None else 7
        enhancement = enhancement if enhancement else "Default"

        request_url = (
            f"{retrosynthesis_base_url}/retrosynthesis_api/"
            f"?key={retrosynthesis_api_key}"
            f"&smiles={validated_smiles}"
            f"&enhancement={enhancement}"
            f"&iterations={iterations}"
            f"&time_limit={time_limit}"
            f"&max_depth={max_depth}"
        )

        # Call the API to start the retrosynthesis job
        unique_identifier = retrosynthesis_process_wrapper(request_url)
        interval = dcc.Interval(id="interval-component", interval=5000, n_intervals=0)

        return (
            "Retrosynthesis process started. Please wait...",
            unique_identifier,
            interval,
            "show",
        )

    @dash_app.callback(
        Output("computed-retrosynthesis-routes", "data"),
        Output("user-message", "children", allow_duplicate=True),
        Output("interval-container", "children", allow_duplicate=True),
        Output("loading-display", "display", allow_duplicate=True),
        Input("interval-component", "n_intervals"),
        State("computed-retrosynthesis-uuid", "data"),
        prevent_initial_call=True,
    )
    def check_retrosynthesis_status(
        n_intervals: int, task_id: str
    ) -> Tuple[Optional[dict], str, Optional[None], str]:
        """
        Check the status of the retrosynthesis process and remove interval element and return results if completed.

        Args:
            n_intervals - calls the function at the interval specified in the element.
            task_id - the UUID of the retrosynthesis task

        Returns:
            retrosynthesis routes as a dict if completed
            a message to give the user feedback
            children for the dcc.Interval container element
            string to set the display setting for the loading wheel.
        """
        # Poll retrosynthesis API for job results
        request_url = (
            f"{retrosynthesis_base_url}/results/{task_id}"
            f"?key={retrosynthesis_api_key}"
        )
        retro_api_status, solved_routes, raw_routes = (
            retrosynthesis_api.retrosynthesis_results_poll(request_url)
        )

        if retro_api_status == "running":
            return dash.no_update, "Processing...", dash.no_update, dash.no_update

        if retro_api_status == "error":
            return dash.no_update, f"Retrosynthesis failed.", None, "hide"

        retrosynthesis_output = {"uuid": task_id, "routes": solved_routes}
        return (
            retrosynthesis_output,
            "Interactive display for retrosynthesis completed.",
            None,
            dash.no_update,
        )

    @dash_app.callback(
        Output("loading-display", "display", allow_duplicate=True),
        Output("computed-conditions-data", "data"),
        State("computed-retrosynthesis-uuid", "data"),
        Input("computed-retrosynthesis-routes", "data"),
        prevent_initial_call=True,
    )
    def new_conditions(
        unique_identifier: str, solved_routes: dict
    ) -> Tuple[Optional[str], Optional[dict]]:
        """
        Called upon completion of a new retrosynthesis routes
        Generates conditions for each corresponding forward reaction in the retrosynthetic routes

        Args:
            # Inputs
            solved_routes - the solved retrosynthetic routes
            # States
            unique_identifier - the uuid for the retrosynthesis

        Returns:
            a string to hide the loading screen if there are no routes to find conditions for
            A dbc table containing the conditions data

        Fires on completion of new retrosynthesis routes
        Inputs: Retrosynthesis routes - using smiles of reactants and products
        makes an api call to get the condition data from the reaction smiles.
        Condition data includes: accuracy score, temperature, solvent, reagents, catalysts
        Returns a Bool to quit the loading circle and a dbc table with conditions data
        """
        if solved_routes and solved_routes != {}:
            conditions = conditions_api.get_conditions(solved_routes["routes"])
            conditions_output = {"uuid": unique_identifier, "routes": conditions}
            return dash.no_update, conditions_output
        return "hide", dash.no_update

    @dash_app.callback(
        Output("active-conditions-data", "data"),
        Input("computed-conditions-data", "data"),
        Input("reloaded-conditions-data", "data"),
        Input("user-uploaded-conditions-data", "data"),
        prevent_initial_call=True,
    )
    def determine_active_conditions(
        computed_conditions: dict,
        reloaded_conditions: dict,
        user_uploaded_conditions: dict,
    ) -> dict:
        """
        Called when one of: new conditions data computed, reloading a retrosynthesis, user uploaded a retrosynthesis
        The context is used to select which one of these conditions sources has changed and is therefore active

        Args:
            # Inputs
            computed_conditions - from the conditions API
            reloaded_conditions - from the database
            user_uploaded_conditions - from a file upload

        Returns:
            the active condition set as a dictionary
        """
        condition_set = ctx.triggered[0]["value"]
        if condition_set:
            return condition_set

    @dash.callback(
        Output("computed-sustainability-data", "data"),
        Input("computed-conditions-data", "data"),
        prevent_initial_call=True,
    )
    def sustainability_assessment(all_conditions: dict) -> dict:
        """
        Called when condition predictions have been processed and outputted to html element.
        Gets the sustainability data for a route from the conditions_dict which contains all the necessary data

        Args:
            all_conditions - conditions data in dictionary

        Returns:
            sustainability_for_all_routes - dictionary of sustainability data
        """
        sustainability_for_all_routes = sustainability.AllRouteSustainability(
            all_conditions["routes"]
        ).get()
        return sustainability_for_all_routes

    @dash.callback(
        Output("active-sustainability-data", "data"),
        Input("computed-sustainability-data", "data"),
        Input("reloaded-sustainability-data", "data"),
        Input("user-uploaded-route-sustainability-data", "data"),
        prevent_initial_call=True,
    )
    def determine_active_sustainability_data(
        sustainability_data: dict,
        reloaded_sustainability: dict,
        user_uploaded_sustainability: dict,
        prevent_initial_call=True,
    ) -> dict:
        """

        Fires on a new retrosynthesis when the conditions api has returned data or when reloading a retrosynthesis.
        Data is used from whichever input triggered the callback
        Returns a list of sustainability data
        """
        active_sustainability_data = ctx.triggered[0]["value"]
        if active_sustainability_data:
            return active_sustainability_data
        return dash.no_update

    @dash_app.callback(
        Output("weighted-sustainability-data", "data"),
        State("active-retrosynthesis-uuid", "data"),
        Input("active-sustainability-data", "data"),
        Input(
            component_id={"type": "sustainability-weighting-slider", "property": ALL},
            component_property="value",
        ),
        prevent_initial_call=True,
    )
    def apply_weightings(
        unique_identifier: str,
        sustainability_data: dict,
        sustainability_weightings: List[int],
    ) -> dict:
        """
        Called when user changes the metric weighting slides or when the active sustainability data changes
        Applies weightings to the sustainability data to generate the weighted sustainability data for each step
        and the route as a whole.

        Args:
            # Inputs
            sustainability - dict of sustainability data
            sustainability_weightings - list of the metric weightings taken from the slider
            # States
            unique_identifier - the uuid for the retrosynthesis

        Returns:
            An updated weighted sustainability dictionary.
        """
        # apply weightings to each route to obtain a weighted median for each route as a whole.
        for route_label, route in sustainability_data.items():
            route_median = sustainability.weighted_median_for_route(
                route, sustainability_weightings
            )
            route["route_average"]["weighted_median"] = route_median
        # 2) apply weighting to each step of the route
        for route in sustainability_data.values():
            sustainability.weighted_median_for_each_step(
                route, sustainability_weightings
            )
        weighted_sustainability_output = {
            "uuid": unique_identifier,
            "routes": sustainability_data,
        }
        return weighted_sustainability_output

    @dash.callback(
        Output("active-retrosynthesis-uuid", "data"),
        Input("computed-retrosynthesis-uuid", "data"),
        Input("reloaded-retrosynthesis-uuid", "data"),
        Input("user-uploaded-route-uuid", "data"),
        prevent_initial_call=True,
    )
    def determine_active_uuid(
        computed_uuid: str, reloaded_uuid: str, user_uploaded_uuid: str
    ) -> str:
        """
        Called when one of: a new retrosynthesis is made, a retrosynthesis is reloaded, retrosynthesis file upload
        Determines the active uuid by using the context

        Args:
            # Inputs
            computed_uuid - for a newly computed retrosynthesis from the smiles input field
            reloaded_uuid - from the database
            user_uploaded_uuid - newly computed after user has uploaded a routes/conditions file

        Returns:
            The active uuid string
        """
        active_uuid_data = ctx.triggered[0]["value"]
        if active_uuid_data:
            return active_uuid_data
        return dash.no_update

    @dash_app.callback(
        Output("active-retrosynthesis-routes", "data"),
        Input("computed-retrosynthesis-routes", "data"),
        Input("reloaded-retrosynthesis-routes", "data"),
        Input("user-uploaded-route", "data"),
        prevent_initial_call=True,
    )
    def determine_active_retrosynthesis_routes(
        computed_retrosynthesis_routes: dict,
        reloaded_retrosynthesis_routes: dict,
        user_route: dict,
    ) -> dict:
        """
        Called when one of: a new retrosynthesis is made, a retrosynthesis is reloaded, retrosynthesis file upload
        Determines the active routes using the context

        Args:
            # Inputs
             computed_retrosynthesis_routes - fresh from the retrosynthesis api
             reloaded_retrosynthesis_routes - from the database
             user_route - from the user uploaded file

        Returns:
            The active retrosynthestic routes
        """
        active_retrosynthesis_data = ctx.triggered[0]["value"]
        if active_retrosynthesis_data:
            return active_retrosynthesis_data
        return dash.no_update

    @dash_app.callback(
        Output("routes-dropdown", "value"),
        Output("routes-dropdown", "options"),
        State("active-retrosynthesis-routes", "data"),
        State("active-conditions-data", "data"),
        Input("weighted-sustainability-data", "data"),
        prevent_intial_call=True,
    )
    def populate_routes_dropdown(
        active_retrosynthesis: dict,
        active_conditions: dict,
        active_weighted_sustainability: dict,
    ) -> Tuple[str, List[dict]]:
        """
        Populates the routes dropdown at the top of the page. Checks all unique identifiers in dictionaries match
        this confirms they are from the same retrosynthesis

        Args:
            # Inputs
            active_retrosynthesis - dictionary with retrosynthesis data of retrosynthesis to display
            active_conditions - dictionary with conditions data of retrosynthesis to display
            active_weighted_sustainability - dictionary with sustainability data of retrosynthesis to display

        Returns:
            Initial route dropdown value defaults to "Route 1"
            List of route options coloured by their sustainability
        """

        if (
            active_retrosynthesis
            and active_conditions
            and active_weighted_sustainability
        ):
            unique_identifier_list = [
                x["uuid"]
                for x in [
                    active_retrosynthesis,
                    active_conditions,
                    active_weighted_sustainability,
                ]
            ]
            # all unique identifiers should match - indicating they come from the same retrosynthesis
            if all(
                unique_identifier == unique_identifier_list[0]
                for unique_identifier in unique_identifier_list
            ):
                route_options = dropdowns.routes(
                    active_retrosynthesis, active_weighted_sustainability
                )
                return (
                    "Route 1",
                    route_options,
                )
        return dash.no_update

    @dash_app.callback(
        Output("routes-dropdown", "style"),
        State("routes-dropdown", "options"),
        Input("routes-dropdown", "value"),
        prevent_initial_call=True,
    )
    def update_route_dropdown_background_colour(
        options: List[dict], active_route: str
    ) -> dict:
        """
        Called when the routes dropdown changes
        Updates the background colour to the weighted median sustainability of the active route in the dropdown

        Args:
            # Inputs
            active_route - the name of the active route in the pattern Route 1, Route 2, etc.
            # States
            options - the options in the dropdown including background colour data

        Returns:
            a dict with the background colour reflective of the sustainability of the selected route
        """
        # get background colour
        for option in options:
            if option["value"] == active_route:
                background_colour = option["label"]["props"]["style"][
                    "background-color"
                ]
                return {"background-color": background_colour, "width": "100%"}

    @dash_app.callback(
        Output("loading-display", "display", allow_duplicate=True),
        Output("retrosynthesis-cytoscape", "elements"),
        Output("retrosynthesis-cytoscape", "stylesheet"),
        State("active-retrosynthesis-routes", "data"),
        Input("routes-dropdown", "value"),
        prevent_initial_call=True,
    )
    def display_retrosynthesis(
        active_retrosynthesis: dict, selected_route: str
    ) -> Tuple[str, List[dict], List[dict]]:
        """
        Called when there is a change to the routes dropdown or active routes
        Create the nodes, edges, and stylesheet to generate the interactive cytoscape

        Args:
            # Inputs
            selected_route - e.g., Route 1
            # States
            active_retrosynthesis - the active retrosynthetic routes

        Returns:
            string to control the display of the loading wheel
            style_sheet - a list of styles as dictionaries. Each has a selector and a style
            elements - a list of nodes and edges as dictionaries. node_id is used to identify nodes and connect nodes

        """
        retro_cytoscape = cytoscape.RetrosynthesisCytoscape(
            active_retrosynthesis, selected_route
        )
        elements = retro_cytoscape.make_cytoscape_elements()
        style_sheet = retro_cytoscape.make_cytoscape_stylesheet()
        return (
            "hide",
            elements,
            style_sheet,
        )

    @dash_app.callback(
        Output("route-feedback", "children"),
        State("active-retrosynthesis-routes", "data"),
        Input("weighted-sustainability-data", "data"),
        Input("routes-dropdown", "value"),  # dropdown
        prevent_initial_call=True,
    )
    def generate_route_table(
        retrosynthesis_data: dict, weighted_sustainability: dict, selected_route: str
    ) -> Tuple[dbc.Table, dbc.Table]:
        """
        Called when there is a change to routes dropdown or the active sustainability data
        Generates the two tables in the routes tab. The brief description table
        and the colour-coded step sustainability analysis table.

        Args:
            # Inputs
            weighted_sustainability - the routes sustainability with metric weightings applied
            selected_route - the title of the active route, e.g., 'Route 1'
            # States
            retrosynthesis_data - the active retrosynthetic routes

        Returns:
            Two tables shown on the 'Routes' tab. One with general data and the other sustainability data.
        """
        if retrosynthesis_data and selected_route:
            route_tables = tables.routes(
                retrosynthesis_data["routes"],
                selected_route,
                weighted_sustainability["routes"],
            )
            return route_tables
        return dash.no_update

    """
    Saving Retrosynthesis Results
    """

    @dash_app.callback(
        Output("save-modal", "is_open"),
        Input("open-save-modal", "n_clicks"),
        Input("close-save-modal", "n_clicks"),
        State("save-modal", "is_open"),
    )
    def toggle_modal(n1: int, n2: int, is_open: bool) -> bool:
        """
        Called when user opens or closes 'Save to Workbook' button
        Toggles the modal window open and shut.
        Args:
            # Inputs
            n1 - Calls function if the open button is clicked
            n2 - Calls function if the closed button is clicked
            # States
            is_open - current open status. False means closed before the user click and vice versa for True.

        Returns:
            bool - opposite of current bool.
        """
        if n1 or n2:
            return not is_open
        return is_open

    @dash_app.callback(
        Output("new-reaction-modal", "is_open"),
        Input("open-new-reaction-modal", "n_clicks"),
        Input("close-new-reaction-modal", "n_clicks"),
        State("new-reaction-modal", "is_open"),
    )
    def toggle_new_reaction_modal(n1: int, n2: int, is_open: bool) -> bool:
        """
        Called when user opens or closes 'Export to Sketcher' button on the Reactions tab.
        Toggles the modal window open and shut
        Args:
            # Inputs
            n1 - Calls function if the open button is clicked
            n2 - Calls function if the closed button is clicked
            # States
            is_open - current open status. False means closed before the user click and vice versa for True.

        Returns:
            bool - opposite of current bool.
        """
        if n1 or n2:
            return not is_open
        return is_open

    @dash_app.callback(
        Output("search-config-collapse", "is_open"),
        Input("search-config-toggle", "n_clicks"),
        State("search-config-collapse", "is_open"),
        prevent_initial_call=True,
    )
    def toggle_search_configuration(n_clicks, is_open):
        """
        Toggles the search configuration panel when the settings button (gear icon) is clicked.

        Args:
            n_clicks (int): The number of times the toggle button is clicked.
            is_open (bool): The current state of the search configuration panel.

        Returns:
            bool: The opposite of the current state (i.e., toggles open/close).
        """
        return not is_open

    @dash_app.callback(
        Output("modal-save-message", "children"),
        Output("new-retrosynthesis-saved-flag", "data"),
        State("save-modal-name", "value"),
        State("active-retrosynthesis-routes", "data"),
        State("active-conditions-data", "data"),
        State("active-sustainability-data", "data"),
        State("save-modal-workbook-dropdown", "value"),
        State("new-retrosynthesis-saved-flag", "data"),
        State("save-functionality-status", "data"),
        State("active-retrosynthesis-uuid", "data"),
        Input("save-modal-save-button", "n_clicks"),
        prevent_initial_call=True,
    )
    def save_retrosynthesis(
        name: str,
        solved_routes: dict,
        conditions: dict,
        sustainability_data: dict,
        workbook_id: int,
        new_retrosynthesis_saved_tracker: int,
        functionality_status: Union[Literal["disabled"], Literal["enabled"]],
        retrosynthesis_uuid: str,
        click_save: int,
    ) -> [str, int]:
        """
        Called when the user clicks 'Save' to save a retrosynthesis
        Validates and saves the current retrosynthesis to the database.
        Args:
            # Inputs:
            click_save - clicking the save button increments the interger, calling this function
            # States:
            name: name of saved retrosynthesis
            solved_routes: retrosynthetic route data dict
            conditions: condition data dict
            sustainability_data: sustainability data dict
            workbook_id: workbook database primary key id
            new_retrosynthesis_saved_tracker: tracker to update saved retrosynthesis list upon change
            functionality_status: values of either 'enabled' or 'disabled' to determine if saving is active

        Returns:
            User message
            tracker int to indicate if changes to the saved reaction list are needed upon successful save.
        """
        if utils.functionality_disabled_check(functionality_status):
            return dash.no_update
        (
            user_message,
            retrosynthesis_saved_tracker,
        ) = saved_retrosyntheses.SaveRetrosynthesis(
            name,
            solved_routes,
            conditions,
            sustainability_data,
            workbook_id,
            new_retrosynthesis_saved_tracker,
            retrosynthesis_uuid,
        ).save_process()
        return user_message, retrosynthesis_saved_tracker

    @dash_app.callback(
        Output("saved-results-list", "children"),
        Input("reload-workbook-dropdown", "value"),
        Input("new-retrosynthesis-saved-flag", "data"),
        State("save-functionality-status", "data"),
        prevent_initial_call=True,
    )
    def show_retrosynthesis_list(
        selected_workbook_id: int,
        new_retrosynthesis_saved: int,
        functionality_status: Union[Literal["disabled"], Literal["enabled"]],
    ):
        """
        Called when there is a change to the workbook dropdown or a new retrosynthesis is saved

        Uses the workbook ID to display an HTML card for each retrosynthesis belonging to that workbook
        in the save list for reload.

        Args:
            # Inputs
            new_retrosynthesis_saved - incremented when a new retrosynthesis is saved and causes function to be called
            # States
            selected_workbook_id - Database ID of the selected workbook in the dropdown
            functionality_status - 'enabled' or 'disabled' to allow/disallow saving related methods

        Returns:
            card_group - A list of retrosynthesis as Dash Bootstrap component HTML Cards.
        """
        if utils.functionality_disabled_check(functionality_status):
            return dash.no_update
        card_group = saved_retrosyntheses.make_retrosynthesis_card_list(
            selected_workbook_id
        )
        return card_group

    @dash_app.callback(
        Output("reloaded-retrosynthesis-routes", "data"),
        Output("reloaded-conditions-data", "data"),
        Output("reloaded-sustainability-data", "data"),
        Output("reloaded-retrosynthesis-uuid", "data"),
        State(
            component_id={"type": "retrosynthesis-reload", "index": ALL},
            component_property="value",
        ),
        Input(
            component_id={"type": "retrosynthesis-reload", "index": ALL},
            component_property="n_clicks",
        ),
        State("save-functionality-status", "data"),
        prevent_initial_call=True,
    )
    def reload_retrosynthesis(
        reload_id_values: List[int],
        reload_button_clicks: List[int],
        functionality_status: Union[Literal["disabled"], Literal["enabled"]],
    ) -> Tuple[dict, dict, dict, str]:
        """
        Called when a 'reload' button is clicked

        Args:
            # Inputs
            reload_button_clicks - The reloaded button which is clicked has a 1 in that index otherwise 0 for unclicked.
            # States
            reload_id_values - The retrosynthesis IDs which could be reloaded
            functionality_status - 'enabled' or 'disabled'

        Returns:
            A dictionary for each the retrosynthesis, conditions, sustainability, and the retrosynthesis uuid
        """
        if utils.functionality_disabled_check(functionality_status):
            return dash.no_update
        # find the element that was clicked - if one was clicked

        if saved_retrosyntheses.assert_button_clicked(reload_button_clicks):
            retrosynthesis_to_reload_id = (
                saved_retrosyntheses.get_retrosynthesis_to_reload_id(reload_id_values)
            )
            (
                retrosynthesis_data,
                condition_data,
                sustainability_data,
                retrosynthesis_uuid,
            ) = saved_retrosyntheses.get_reloaded_retrosynthesis(
                retrosynthesis_to_reload_id
            )
            return (
                retrosynthesis_data,
                condition_data,
                sustainability_data,
                retrosynthesis_uuid,
            )
        return dash.no_update

    """
    Compound sidebar
    """

    @dash_app.callback(
        Output("compound-feedback", "children"),
        Output("tapped-compound-image", "src"),
        Input("retrosynthesis-cytoscape", "tapNodeData"),
        prevent_initial_call=True,
    )
    def display_compound_node_data(
        tapped_node: dict,
    ) -> Tuple[Union[dbc.Table, str], str]:
        """
        Called when user taps a compound node in the cytoscape interface.
        Uses the compound smiles to make an image and find the compound in the database if it is present

        Args:
            # Inputs
            tapped_node - corresponds to a molecule with the SMILES inside the dictionary

        Returns:
            either the compound data table or a string stating compound is not in the database
            a svg image of the compound
        """
        if tapped_node == ["smiles"] or tapped_node is None:
            return dash.no_update
        smiles = tapped_node["smiles"]
        img_data = utils.rdkit_smiles_to_image(smiles)
        compound_feedback = tables.compound(smiles)
        return compound_feedback, img_data

    """
    Conditions
    """

    @dash_app.callback(
        Output("reaction-conditions-list", "data"),
        Output("reaction-sustainability-list", "data"),
        Output("conditions-dropdown", "value"),  # dropdown
        Output("conditions-dropdown", "options"),
        State("routes-dropdown", "value"),
        State("active-conditions-data", "data"),
        Input("retrosynthesis-cytoscape", "tapNodeData"),
        Input("weighted-sustainability-data", "data"),
        prevent_initial_call=True,
    )
    def fill_conditions_dropdown(
        route_label: str,
        conditions_data: Dict,
        tapped_node: Dict,
        weighted_sustainability_data: Dict,
    ) -> Tuple[List[Dict], List[Dict], str, List[Dict]]:
        """
        Called when a user clicks on a compound node
        Finds the data for all the condition sets (up to 10) for a forward reaction and makes the dropdown for this.

        Args:
            # Inputs
            tapped_node - The node the user has clicked on. product SMILES used to look up current reaction
            # States
            route - The current route label - needed to look up current route
            conditions_data - The reaction conditions are extracted from this dictionary

        Returns:
            rxn_conditions - conditions for the current reaction
            rxn_sustainability - weighted sustainability for the current reaction
            'Condition Set 1' as the default active condition set
            dropdown_options - to populate the condition set dropdown.

        """
        if (
            tapped_node
            and tapped_node.get("reaction_smiles")
            and tapped_node["uuid"] == conditions_data["uuid"]
        ):
            (
                rxn_conditions,
                rxn_sustainability,
                dropdown_options,
            ) = dropdowns.make_conditions_dropdown(
                route_label,
                conditions_data["routes"],
                weighted_sustainability_data["routes"],
                tapped_node,
            )
            return (
                rxn_conditions,
                rxn_sustainability,
                "Condition Set 1",
                dropdown_options,
            )
        # a terminal node has no reaction SMILES so return empty values
        else:
            return (
                [{}],
                [{}],
                "Terminal node.",
                dropdowns.make_terminal_node_dropdown(),
            )

    @dash_app.callback(
        Output("conditions-dropdown", "style"),
        State("conditions-dropdown", "options"),
        Input("conditions-dropdown", "value"),
        prevent_initial_call=True,
    )
    def update_conditions_dropdown_background_colour(
        options: List[dict], active_condition_set: str
    ) -> dict:
        """
        Called when the conditions dropdown changes
        Updates the background colour to the weighted median sustainability of the active condition set in the dropdown

        Args:
            # Inputs
            active_condition set - the name of the active condition set in the pattern Condition Set 1, Condition Set 2, etc.
            # States
            options - the options in the dropdown including background colour data

        Returns:
            a dict with the background colour reflective of the sustainability of the selected condition set
        """

        if active_condition_set != "Terminal node":
            # get background colour
            for option in options:
                if option["value"] == active_condition_set:
                    background_colour = option["label"]["props"]["style"][
                        "background-color"
                    ]
                    return {
                        "background-color": background_colour,
                        "width": "100%",
                        "margin-bottom": "1rem",
                    }
        # if a terminal node return empty values
        else:
            return {
                "background-color": "#FFFFFF",
                "width": "100%",
                "margin-bottom": "1rem",
            }

    @dash_app.callback(
        Output("reaction-conditions", "children"),
        # Output('reaction-sustainability', 'children'),
        State("reaction-conditions-list", "data"),
        State("reaction-sustainability-list", "data"),
        Input("conditions-dropdown", "value"),
        prevent_initial_call=True,
    )
    def generate_reaction_table(
        conditions_options: Dict[str, dict],
        sustainability_options: Dict[str, dict],
        conditions_dropdown_value: str,
    ) -> Union[dbc.Table, str]:
        """
        Called when the conditions dropdown changes or clicks a chemical node to show details of the predicted reaction
        Generates the reaction table in the reactions tab with the details predicted to perform the forward reaction
        colour coded by their sustainability.

        Args:
            # Inputs
            conditions_dropdown_value - The label/value of the active condition set in format: 'Condition Set 1'
            # States
            conditions_options - Dict with the conditions for the forward reaction to be shown in the table
            sustainability_options - Dict with the sustainability for the conditions - colours the rows in the table.

        Returns:
            Either the colour-coded conditions table for the forward reaction or a string explaining terminal node
            has no reaction.
        """
        if conditions_dropdown_value and conditions_dropdown_value != "Terminal node.":
            conditions_table = tables.reaction(
                conditions_options, sustainability_options, conditions_dropdown_value
            )
            return conditions_table
        elif conditions_dropdown_value:
            return "Terminal node has no reaction"
        return dash.no_update

    """Route data"""

    @dash_app.callback(
        Output("reaction-smiles", "data"),
        Output("tapped-reaction-image", "src"),
        Output("reaction-class", "children"),
        State("window-width", "data"),
        Input("retrosynthesis-cytoscape", "tapNodeData"),
        prevent_initial_call=True,
    )
    def display_reaction(window_width: int, tapped_node: dict) -> Tuple[str, str, str]:
        """
        Called when user clicks on a compound node
        Uses the reaction string to make a png image with the reaction class above the image in the Reactions tab.

        Args:
            # Inputs
            tapped_node - dictionary contains reaction_smiles and reaction_class of the active reaction

        Returns:
            reaction_smiles of the active reaction
            img_data for the current reaction as a png string
            reaction_class of the active reaction.
        """
        reaction_class = tapped_node["label"]
        reaction_smiles = tapped_node["reaction_smiles"]
        if reaction_smiles:
            img_data = utils.reaction_smiles_to_image(window_width, reaction_smiles)
        else:
            img_data = utils.rdkit_smiles_to_image(
                tapped_node["smiles"]
            )  # rdkit method best for single singles
        return reaction_smiles, img_data, reaction_class

    @dash_app.callback(
        Output("new-reaction-id", "value"),
        State("save-functionality-status", "data"),
        Input("new-reaction-workbook-dropdown", "value"),
        prevent_initial_call=True,
    )
    def update_new_reaction_id(
        functionality_status: Union[Literal["disabled"], Literal["enabled"]],
        workbook_id: int,
    ) -> str:
        """
        Called when changing the workbook dropdown and finds the next reaction ID. This is needed to make a new reaction

        Args:
            # Inputs
            workbook_id - the ID of the selected workbook in the dropdown
            # States
            functionality_status - saving methods are 'enabled' or 'disabled'

        Returns:
            the next auto-incremented reaction_id value as a string. e.g., WB1-001

        """
        if utils.functionality_disabled_check(functionality_status) or not workbook_id:
            return dash.no_update
        return services.reaction.get_next_reaction_id_for_workbook(workbook_id)

    @dash_app.callback(
        Output("new-reaction-url", "data"),
        Output("modal-new-reaction-message", "children"),
        State("url", "pathname"),
        State("new-reaction-workbook-dropdown", "value"),
        State("new-reaction-name", "value"),
        State("new-reaction-id", "value"),
        State("reaction-smiles", "data"),
        State("save-functionality-status", "data"),
        Input("new-reaction-data-submit", "n_clicks"),
        prevent_initial_call=True,
    )
    def new_reaction(
        current_url: str,
        workbook_id: int,
        reaction_name: str,
        reaction_id: str,
        reaction_smiles: str,
        functionality_status: Union[Literal["disabled"], Literal["enabled"]],
        n_clicks: int,
    ) -> Tuple[str, str]:
        """
        Called when user clicks new reaction button - to export reaction from retrosynthesis to new ELN entry
        Performs validation and then if successful saves Reaction to database and opens in new tab

        Args:
            # Inputs
            n_clicks - integer changes when user clicks the new reaction button
            # States
            current_url - to get the base url
            workbook_id - ID used to save the new reaction
            reaction_name - save to database under this name
            reaction_id - saved to database under this id
            reaction_smiles - saved to database
            functionality_status - 'enabled' or 'disabled'

        Returns:
            URL of the new reaction and opens in a new tab
            A feedback message to user in case of failure

        """
        if utils.functionality_disabled_check(functionality_status):
            return dash.no_update
        if not reaction_smiles:
            return "", "Cannot make a new reaction from a terminal node."
        # get workbook and make a reaction from it
        workbook_object = services.workbook.get(workbook_id)
        result = saved_retrosyntheses.save_new_reaction_from_retrosynthesis(
            workbook_object, reaction_name, reaction_id, reaction_smiles
        )
        # if this is successful make the url and redirect the user
        if result == "New reaction made":
            workgroup_name = workbook_object.WorkGroup.name
            workbook_name = workbook_object.name
            base_url = current_url.split("retrosynthesis")[0]
            new_url = quote(
                f"{base_url}sketcher/{workgroup_name}/{workbook_name}/{reaction_id}/no"
            )
            return new_url, "New reaction made!"
        return "", result

    @dash_app.callback(
        Output("url", "pathname"),
        State("save-functionality-status", "data"),
        Input("new-reaction-url", "data"),
        prevent_initial_call=True,
    )
    def go_to_new_reaction(
        functionality_status: Union[Literal["disabled"], Literal["enabled"]],
        new_url: str,
    ) -> str:
        """
        Called when a new reaction is successfully made
        Opens a new reaction in a new tab. Checks url is present first.

        Args:
            functionality_status: 'enabled' or 'disabled'
            new_url: the url which is opened in a new tab for the new ELN reaction entry.
        """
        if utils.functionality_disabled_check(functionality_status):
            return dash.no_update
        if new_url:
            return new_url
        return dash.no_update

    """
    Callbacks for user inputted routes
    """

    @dash_app.callback(
        Output("user-uploaded-route", "data"),
        Output("user-uploaded-conditions-data", "data"),
        Output("user-uploaded-route-uuid", "data"),
        Output("uploaded-route-message", "children"),
        Input("upload-route-button", "contents"),
        State("upload-route-button", "filename"),
        prevent_initial_call=True,
    )
    def upload_user_route_file(
        contents: str, filename: str
    ) -> Tuple[dict, dict, str, str]:
        """
        Called when a user clicks 'Upload Route'
        Uploads the route from the user selected file and shows in the cytoscape

        Args:
            contents - the file contents must be csv, xls, or ods.
            filename - the name of the file the user has uploaded

        Returns:
            processed_route - the dict with route data from the uploaded file
            processed_conditions - the dict with condition data from the uploaded file
            uuid for the uploaded retrosynthesis.
            user feedback - string indicating upload status
        """

        try:
            if contents is not None:
                (
                    processed_route,
                    processed_conditions,
                ) = user_uploaded_routes.read_user_route_file(contents, filename)

                return (
                    processed_route,
                    processed_conditions,
                    processed_route["uuid"],
                    "File successfully uploaded",
                )
        except Exception as e:
            print(e)
            return (
                {},
                {},
                "",
                "Error with the file upload.This may occur if the file extension was changed or an incorrect format "
                "was used.\nPlease refer to the example file for the correct format and upload as .csv, .xslx or .ods.",
            )
        return dash.no_update

    @dash_app.callback(
        Output("user-uploaded-route-sustainability-data", "data"),
        Input("user-uploaded-conditions-data", "data"),
        prevent_initial_call=True,
    )
    def process_user_route_sustainability(conditions: dict) -> dict:
        """
        Get the sustainability data for a user uploaded route after the UUID is available

        Args:
            conditions - the dict with condition data from the uploaded file

        Returns:
            sustainability - the dict with sustainability data from the uploaded file
        """
        sustainability_data = sustainability.AllRouteSustainability(
            conditions["routes"]
        ).get()
        return sustainability_data

    @dash_app.callback(
        Output("example-route-file-download", "data"),
        Input("btn-example-route-file", "n_clicks"),
        prevent_initial_call=True,
    )
    def download_example_file(n_clicks):
        """Downloads the example route file when user"""
        return dcc.send_file(
            os.path.join(
                current_app.config["APP_DIRECTORY"],
                "sources",
                "static",
                "retrosynthesis-route-example.csv",
            )
        )

    dash_app.clientside_callback(
        """
        function(href) {
            var w = window.innerWidth;

            return w;
        }
        """,
        Output("window-width", "data"),
        Input("url", "href"),
    )
    return dash_app.server
