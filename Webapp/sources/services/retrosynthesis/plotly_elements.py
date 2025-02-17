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
            children=html.Div(id="loading-output-1"),
            type="circle",
            fullscreen=True,
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
Enhancement Dropdown with Tooltips
"""
enhancement_dropdown = html.Div(
    className="input-group",
    children=[
        # Dropdown
        html.Div(
            dcc.Dropdown(
                id="enhancement-dropdown",
                options=[
                    {"label": "Default", "value": "Default"},
                    {"label": "Enhanced Search", "value": "eUCT"},
                    {"label": "Enhanced Speed", "value": "dUCT-v1"},
                    {"label": "Enhanced Solve Rate", "value": "dUCT-v2"},
                ],
                placeholder="Select Enhancement Type",
                style={"width": "100%", "border": "1px solid #9FA6B2"},
            ),
            style={"flex-grow": "1"},
        ),
        html.Div(
            html.Span(
                "i",
                id="enhancement-tooltip",
                style={
                    "display": "inline-flex",
                    "align-items": "center",
                    "justify-content": "center",
                    "width": "20px",
                    "height": "20px",
                    "borderRadius": "50%",
                    "backgroundColor": "#007BFF",
                    "color": "white",
                    "fontWeight": "bold",
                    "cursor": "pointer",
                    "margin-left": "8px",
                },
            ),
            style={"display": "flex", "align-items": "center"},
        ),
        # Tooltip Content
        dbc.Tooltip(
            "Options Explained:\n"
            "- Default: Standard configuration for retrosynthesis.\n"
            "- Enhanced Search: Small speed and solve rate improvements\n"
            "- Enhanced Speed: Prioritises speed of results.\n"
            "- Enhanced Solve Rate: Most efficient at solving difficult molecules",
            target="enhancement-tooltip",
            placement="right",
            style={"white-space": "pre-line", "max-width": "250px"},
        ),
    ],
    style={"display": "flex", "align-items": "center", "gap": "5px"},
)

"""
MCTS Search configurations
"""

search_configurations = dbc.Collapse(
    [
        html.Label("Number of Iterations:"),
        dcc.Input(
            id="iterations-input", type="number", value=100, min=20, max=5000, step=1
        ),
        html.Label(" Time Limit (seconds):"),
        dcc.Input(
            id="time-limit-input", type="number", value=60, min=30, max=1000, step=1
        ),
        html.Label(" Max Depth:"),
        dcc.Input(id="max-depth-input", type="number", value=7, min=1, max=25, step=1),
    ],
    id="search-config-collapse",
    is_open=False,
    style={"margin-bottom": "15px", "border": "1px solid #9FA6B2", "padding": "10px"},
)

search_toggle_button = html.Div(
    [
        dbc.Button(
            html.I(className="fa fa-sliders"),
            id="search-config-toggle",
            className="btn btn-outline-secondary btn-sm",
            style={
                "width": "32px",
                "height": "32px",
                "border-radius": "50%",
                "display": "flex",
                "align-items": "center",
                "justify-content": "center",
                "padding": "0",
                "margin-left": "10px",
                "font-size": "16px",
            },
        ),
        dbc.Tooltip(
            "Search Configuration",
            target="search-config-toggle",
            placement="bottom",
            delay={"show": 500, "hide": 100},
        ),
    ]
)

"""
Smiles field + Retrosynthesis button
"""
smiles_field_and_retrosynthesis_button = dbc.Row(
    className="gx-2",
    children=[
        # Column 1: New Retrosynthesis Button and SMILES Input
        dbc.Col(
            width=3,
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
            ],
        ),
        # Column 2: Route Dropdown
        dbc.Col(
            width=3,
            children=[
                inputs,
                html.Div(id="route-score"),
            ],
        ),
        # Column 3: Enhancement Dropdown
        dbc.Col(
            width=3,
            children=[enhancement_dropdown],
        ),
        # Column 4: Save to Workbook Button
        dbc.Col(
            width=3,
            children=[
                html.Button(
                    "Save To Workbook",
                    id="open-save-modal",
                    disabled=False,
                    n_clicks=0,
                    className="btn btn-outline-primary m-1",
                )
            ],
        ),
        dbc.Col(
            width=12,  # Takes full width below the above fields
            children=[
                html.Div(
                    style={"display": "flex", "align-items": "center"},
                    children=[
                        search_toggle_button,
                        search_configurations,
                    ],
                ),
            ],
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


"""Tabs layout and about retrosynthesis button"""
tabs = html.Div(
    children=[
        html.Div(
            html.A(
                "About Retrosynthesis ❯",
                href="/retrosynthesis_about",
                className="btn btn-primary",
                style={
                    "border-radius": "7px",
                    "padding": "10px 14px",
                    "border": "none",
                    "display": "inline-block",
                    "text-decoration": "none",
                },
            )
        ),
        dbc.Tabs(
            [
                dbc.Tab(compound_sidebar, label="Compounds", id="compound-tab"),
                dbc.Tab(reaction_sidebar, label="Reactions", id="reaction-tab"),
                dbc.Tab(route_sidebar, label="Routes", id="route-tab"),
                dbc.Tab(
                    saved_results_sidebar, label="Saved Results", id="saved-results-tab"
                ),
            ]
        ),
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
