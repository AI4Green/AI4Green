import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
from dash import dcc, html

from .style_sheets import WebElementsStyles as webStyle

cyto.load_extra_layouts()

"""
Canvas
"""
retro_tree = cyto.Cytoscape(
    id="retrosynthesis-cytoscape",
    style={"width": "100%", "height": "60rem"},
    layout={
        "name": "dagre",
        "roots": '[id = "node-0"]',
        "fit": True,
        "spacingFactor": 1,
    },
    zoomingEnabled=True,
    responsive=True,
    minZoom=0.5,
    maxZoom=3,
    panningEnabled=True,
    zoom=0.5,
)

"""
Loading circle
"""
loading_circle = html.Div(
    id="loading-div",
    children=[
        dcc.Loading(
            id="loading-display",
            # className="",
            children=html.Div(id="loading-output-1"),
            type="spinner",
            fullscreen=True,
            # loading_state={
            #     "component_name": "retrosynthesis-tree",
            #     # "is_loading": False,
            #     "prop_name": "retro",
            # },
        ),
    ],
)

"""
Dropdowns
"""
inputs = html.Div(
    className="input-group",
    children=[
        html.Div(
            className="input-group-prepend",
            style={"width": "100%"},
            children=[
                dcc.Dropdown(
                    id="routes-dropdown",
                    placeholder="Route (Number of Steps)",
                    style={"width": "100%", "border": "1px solid #9FA6B2"},
                ),
                html.Div(id="routes_output"),
            ],
        )
    ],
    style={"width": "auto"},
)
"""
Smiles field + Retrosynthesis button
"""
smiles_field_and_retrosynthesis_button = dbc.Row(
    className="gx-2",
    children=[
        dbc.Col(
            children=[
                html.Div(
                    className="input-group",
                    children=[
                        html.Div(
                            className="input-group-prepend",
                            children=[
                                html.Button(
                                    "New Retrosynthesis",
                                    id="btn-retrosynthesis",
                                    n_clicks=0,
                                    className="btn btn-outline-primary",
                                )
                            ],
                        ),
                        dcc.Input(
                            id="smiles-input",
                            type="text",
                            className="form-control",
                            placeholder="SMILES",
                        ),
                    ],
                ),
            ]
        ),
        dbc.Col(children=[inputs, html.Div(id="route-score")]),
        dbc.Col(
            children=[
                html.Button(
                    "Save To Workbook",
                    id="open-save-modal",
                    disabled=False,
                    n_clicks=0,
                    className="btn btn-outline-primary m-1",
                )
            ]
        ),
    ],
)

"""Header and inputs"""
header_and_inputs = html.Div(
    children=[
        html.H2("Retrosynthesis"),
        smiles_field_and_retrosynthesis_button,
        # user_route,
        html.H6(
            "Interactive Retrosynthesis Display", id="user-message", className="mt-3"
        ),
        html.Div([loading_circle]),
    ]
)

"""Hidden inputs for data storage"""
retrosynthesis_data = [
    dcc.Store(id="validated-smiles", data="", storage_type="memory"),
    dcc.Store(id="computed-retrosynthesis-routes", storage_type="memory"),
    dcc.Store(id="reloaded-retrosynthesis-routes", data="", storage_type="memory"),
    dcc.Store(id="active-retrosynthesis-routes", data="", storage_type="memory"),
    dcc.Store(id="new-retrosynthesis-saved-flag", data=0, storage_type="memory"),
    dcc.Store(id="computed-retrosynthesis-uuid", data="", storage_type="memory"),
    dcc.Store(id="reloaded-retrosynthesis-uuid", data="", storage_type="memory"),
    dcc.Store(id="active-retrosynthesis-uuid", data="", storage_type="memory"),
    dcc.Store(id="user-uploaded-route-uuid", storage_type="memory"),
    dcc.Store(id="user-uploaded-route", data="", storage_type="memory"),
]

conditions_data = [
    dcc.Store(id="computed-conditions-data", storage_type="memory"),
    dcc.Store(id="reloaded-conditions-data", storage_type="memory"),
    dcc.Store(id="active-conditions-data", storage_type="memory"),
    dcc.Store(id="user-uploaded-conditions-data", storage_type="memory"),
]

sustainability_data = [
    dcc.Store(id="computed-sustainability-data", storage_type="memory"),
    dcc.Store(id="reloaded-sustainability-data", storage_type="memory"),
    dcc.Store(id="active-sustainability-data", storage_type="memory"),
    dcc.Store(id="weighted-sustainability-data", storage_type="memory"),
    dcc.Store(id="user-uploaded-route-sustainability-data", storage_type="memory"),
]

other_data_stores = [
    dcc.Location(id="url"),
    dcc.Store(id="window-width", data=0, storage_type="memory"),
    dcc.Store(id="test-temp", data=0, storage_type="memory"),
    dcc.Store(id="new-reaction-url", data="", storage_type="memory"),
    dcc.Store(id="new-reaction-success", data=False, storage_type="memory"),
    dcc.Store(id="save-functionality-status", data=""),
]

all_data_stores = (
    retrosynthesis_data + conditions_data + sustainability_data + other_data_stores
)

data_storage_elements = html.Div(children=all_data_stores)

"""Compound sidebar"""
compound_sidebar = html.Div(
    [
        html.H3("Compounds"),
        html.Hr(),
        html.Img(id="tapped-compound-image"),
        html.Div(id="compound-feedback"),
    ],
    id="compound-sidebar",
    style=webStyle.SIDEBAR_STYLE,
)

"""Reaction  sidebar"""
reaction_sidebar = html.Div(
    [
        html.H3("Reactions"),
        html.Hr(),
        html.Div(id="reaction-class"),
        html.Img(id="tapped-reaction-image"),
        dcc.Store(id="reaction-smiles", data="", storage_type="memory"),
        html.Div(id="reaction-feedback"),
        dcc.Store(id="reaction-conditions-list", data="", storage_type="memory"),
        html.Div(id="interval-container"),
        dcc.Store(id="reaction-sustainability-list", data="", storage_type="memory"),
        dcc.Dropdown(id="conditions-dropdown"),
        html.Div(id="reaction-conditions"),
        html.Div(id="reaction-sustainability"),
        html.Button(
            "Export to Sketcher",
            id="open-new-reaction-modal",
            className="btn btn-outline-secondary pull-left",
            n_clicks=0,
        ),
    ],
    id="reaction-sidebar",
    style=webStyle.SIDEBAR_STYLE,
)

"""Route sidebar"""
slider_marks = {5: "Default"}


def slider_maker(sustainability_property):
    return html.Div(
        id=f"{sustainability_property}-slider-div",
        style={"display": "grid", "grid-template-columns": "30% 70%"},
        children=[
            html.Label(f'{sustainability_property.replace("-", " ").title()}'),
            dcc.Slider(
                0,
                10,
                1,
                value=5,
                marks=slider_marks,
                tooltip={"placement": "bottom", "always_visible": False},
                id={
                    "type": "sustainability-weighting-slider",
                    "property": sustainability_property,
                },
            ),
        ],
    )


route_sidebar = html.Div(
    [
        html.H3("Routes"),
        html.Hr(),
        html.Div(
            [
                dcc.Upload(
                    "Upload Route",
                    className="btn btn-outline-secondary m-1",
                    id="upload-route-button",
                ),
                html.Button(
                    "Example Route File",
                    id="btn-example-route-file",
                    className="btn btn-outline-secondary m-1",
                ),
                dcc.Download(id="example-route-file-download"),
            ],
            style={"display": "flex", "justify-content": "space-between"},
        ),
        html.P(id="uploaded-route-message"),
        html.Hr(),
        html.Div(id="route-feedback"),
        html.Div(id="sustainability-feedback"),
        html.Div(
            id="sustainability-sliders-div",
            className="pb-5",
            children=[
                html.Hr(),
                html.H6("Adjust Sustainability Metric Weightings"),
                html.Hr(),
                slider_maker("solvent"),
                slider_maker("temperature"),
                slider_maker("stoichiometry/catalyst"),
                slider_maker("elements"),
                slider_maker("atom-economy"),
                slider_maker("safety"),
            ],
        ),
    ],
    id="route-sidebar",
    style=webStyle.SIDEBAR_STYLE,
)

"""Saved results sidebar"""
saved_results_sidebar = html.Div(
    [
        html.H3("Saved Results"),
        dcc.Dropdown(
            id="reload-workbook-dropdown", className="pb-3", placeholder="workbook"
        ),
        html.Hr(),
        html.Label("Saved results list"),
        html.Hr(),
        html.Div(id="saved-results-list", className="row d-flex justify-content-left"),
        html.Data(id="are-images-shown", value=True),
    ],
    id="saved-results-sidebar",
    style=webStyle.SIDEBAR_STYLE,
)

tabs = dbc.Tabs(
    [
        dbc.Tab(compound_sidebar, label="Compounds", id="compound-tab"),
        dbc.Tab(reaction_sidebar, label="Reactions", id="reaction-tab"),
        dbc.Tab(route_sidebar, label="Routes", id="route-tab"),
        dbc.Tab(saved_results_sidebar, label="Saved Results", id="saved-results-tab"),
    ]
)
"""Save results modal window"""

save_modal = html.Div(
    [
        dbc.Modal(
            [
                dbc.ModalHeader(
                    dbc.ModalTitle("Save Retrosynthesis"), close_button=False
                ),
                dbc.ModalBody(
                    children=[
                        dbc.Label("Workbook", html_for="save-modal-workbook-dropdown"),
                        dcc.Dropdown(
                            id="save-modal-workbook-dropdown",
                            className="pb-3",
                            placeholder="workbook",
                        ),
                        dbc.Label("Name", html_for="save-modal-name"),
                        dbc.Input(
                            id="save-modal-name", className="pb-3", placeholder="name"
                        ),
                        html.Div(id="modal-save-message")
                        # html.Div(
                        #
                        # )
                    ]
                ),
                dbc.ModalFooter(
                    children=[
                        dbc.Button(
                            "Close",
                            id="close-save-modal",
                            className="btn-danger",
                            n_clicks=0,
                        ),
                        dbc.Button(
                            "Save",
                            id="save-modal-save-button",
                            className="btn btn-success",
                            n_clicks=0,
                        ),
                    ]
                ),
            ],
            id="save-modal",
            is_open=False,
        ),
    ]
)

new_reaction_modal = html.Div(
    [
        dbc.Modal(
            [
                dbc.ModalHeader(dbc.ModalTitle("New Reaction"), close_button=False),
                dbc.ModalBody(
                    children=[
                        html.Form(
                            [
                                dbc.Label(
                                    "Workbook",
                                    html_for="new-reaction-workbook-dropdown",
                                ),
                                dcc.Dropdown(
                                    id="new-reaction-workbook-dropdown",
                                    className="pb-3",
                                    placeholder="workbook",
                                ),
                                dbc.Label("Reaction ID", html_for="new-reaction-id"),
                                dbc.Input(
                                    id="new-reaction-id",
                                    className="pb-3",
                                    disabled=True,
                                    type="text",
                                    size="sm",
                                ),
                                dbc.Label(
                                    "Reaction Name", html_for="new-reaction-name"
                                ),
                                dbc.Input(
                                    id="new-reaction-name",
                                    className="pb-3",
                                    type="text",
                                    size="sm",
                                ),
                            ]
                        ),
                        html.Div(id="modal-new-reaction-message"),
                    ]
                ),
                dbc.ModalFooter(
                    children=[
                        dbc.Button(
                            "Cancel",
                            id="close-new-reaction-modal",
                            className="btn-danger",
                            n_clicks=0,
                        ),
                        dbc.Button(
                            "Create",
                            id="new-reaction-data-submit",
                            className="btn-success",
                            n_clicks=0,
                        ),
                    ]
                ),
            ],
            id="new-reaction-modal",
            is_open=False,
        ),
    ]
)

references = html.Div(
    [
        html.Br(),
        html.B("References:"),
        html.P(
            [
                "Retrosynthetic routes predicted with: AiZynthFinder: a fast, robust and flexible open-source software for "
                "retrosynthetic planning. ",
                html.I("J. Cheminform."),
                ", 2020, ",
                html.B("12"),
                ", 70 (2020).",
            ]
        ),
        html.P(
            [
                "Reaction conditions predicted with: Using Machine Learning To Predict Suitable Conditions for Organic "
                "Reactions. ",
                html.I("ACS Cent. Sci."),
                ", 2018, ",
                html.B("4"),
                ", 11, 1465–1476. ",
            ]
        ),
    ]
)
