from typing import Dict, List, Optional, Tuple, Union

import dash_bootstrap_components as dbc
import pandas as pd
from dash import html
from flask_login import current_user
from sources import models, services


def compound(smiles: str) -> Union[dbc.Table, str]:
    """
    Finds the compound in the database via InCHI and makes a table if it is present

    Args:
        smiles - the SMILES of the compound of interest
    Returns:
        Either a compound data table or a string stating compound is not in the database.
    """
    inchi = services.all_compounds.smiles_to_inchi(smiles)
    if inchi is None:
        return "Cannot parse structure to InChI"
    compound_object = services.compound.get_compound_from_inchi(inchi)
    if compound_object is None:
        return f"Compound with SMILES {smiles} is not in the database"
    compound_table = make_compound_table(compound_object)
    return compound_table


def make_compound_table(compound_object: models.Compound) -> dbc.Table:
    """
    We make a compound table for the compound tab.

    Args:
        compound - from the database

    Returns:
        A table with compound name mol. weight, CAS and hazards.
    """
    df = pd.DataFrame(
        {
            "": ["Name", "Molecular Weight", "CAS", "Hazard Codes"],
            "Compound": [
                compound_object.name,
                compound_object.molec_weight,
                compound_object.cas,
                compound_object.hphrase,
            ],
        }
    )
    return dbc.Table.from_dataframe(
        df, bordered=True, hover=True, id="compound-data-table"
    )


def reaction(
    conditions_options: Dict[str, dict],
    sustainability_options: Dict[str, dict],
    condition_set_label: str,
) -> dbc.Table:
    """
    Generates a table for reaction conditions and sustainability.

    Args:
        conditions_options: The conditions for the reaction
        sustainability_options: The sustainability for the reaction conditions - controls background colour of cells
        condition_set_label: The selected condition set

    Returns:
        dbc.Table: Table containing reaction conditions and sustainability information.
    """
    # selected_condition_idx = (
    #         int(condition_set_label[-1]) - 1
    # )  # minus 1 to use zero-based index
    conditions_data = conditions_options[condition_set_label]
    sustainability = sustainability_options[condition_set_label]
    # update both dictionaries in preparation of making the table
    conditions_data = format_conditions(conditions_data)
    get_table_cell_styles(sustainability)
    table_header = [
        html.Thead(
            html.Tr(
                [
                    html.Th("", style={"border-color": "black"}),
                    html.Th(
                        "Conditions & Sustainability", style={"border-color": "black"}
                    ),
                ]
            )
        )
    ]
    table_body = make_reaction_table_rows(conditions_data, sustainability)
    table = dbc.Table(
        table_header + table_body,
        bordered=True,
        id="reaction-conditions-data-table",
    )
    return table


def make_reaction_table_rows(
    conditions_data: Dict, sustainability: Dict
) -> List[html.Tbody]:
    """
    Generates a table body for conditions and sustainability information.

    Args:
        conditions_data (dict): Data on reaction conditions.
        sustainability (dict): Sustainability information.

    Returns:
        list: List containing table body elements.
    """
    temperature_row = generate_html_row(
        "Temperature (°C)",
        conditions_data["temperature"],
        sustainability["temperature"]["style"],
    )
    solvent_row = generate_html_row(
        "Solvent",
        ", ".join(conditions_data["solvent_names"]),
        sustainability["solvent"]["style"],
    )
    reagents_row = generate_html_row(
        "Reagents", ", ".join(conditions_data["reagents_names"])
    )
    catalyst_row = generate_html_row(
        "Stoichiometry/Catalyst",
        ", ".join(conditions_data["catalyst_names"]),
        sustainability["catalyst"]["style"],
    )
    element_sustainability_row = generate_html_row(
        "Element Sustainability",
        sustainability["element_sustainability"]["years"],
        sustainability["element_sustainability"]["style"],
    )
    atom_economy_row = generate_html_row(
        "Atom Economy",
        f"{sustainability['atom_economy']['atom_economy']}%",
        sustainability["atom_economy"]["style"],
    )
    if len(sustainability["safety"]["hazard_code_list"]) == 0:
        safety_text = "No hazard codes found"
    else:
        safety_text = ", ".join(sustainability["safety"]["hazard_code_list"])
    safety_row = generate_html_row(
        "Safety",
        safety_text,
        sustainability["safety"]["style"],
    )

    table_body = [
        html.Tbody(
            [
                temperature_row,
                solvent_row,
                reagents_row,
                catalyst_row,
                element_sustainability_row,
                atom_economy_row,
                safety_row,
            ]
        )
    ]
    return table_body


def get_table_cell_styles(sustainability: Dict):
    """
    Updates styles for sustainability table cells.

    Args:
        sustainability (dict): Sustainability information
    """
    hazard_colours = current_user.hazard_colors
    flag_score_to_phrase_dict = {
        4: "HighlyHazardous",
        3: "Hazardous",
        2: "Problematic",
        1: "Recommended",
    }
    for key, value in sustainability.items():
        if type(value) is dict:
            score = value["flag"]
            if isinstance(score, float):
                score = round(score, 0)
            chem21_phrase = flag_score_to_phrase_dict[score]
            background_colour = hazard_colours[chem21_phrase]
            text_colour = hazard_colours[f"{chem21_phrase}_text"]
            style_dict = {
                "background-color": background_colour,
                "color": text_colour,
                "border-color": "black",
                "word-wrap": "break-word",
                "white-space": "normal",
            }
            sustainability[key]["style"] = style_dict


def format_conditions(conditions_data: Dict) -> Dict:
    """
    Formats condition data, replacing the null values of "Not Found" and "No Compound" with the SMILEs and "None"
    respectively

    Args:
        conditions_data: The dictionary with the conditions data including reagents, catalyst, amd solvent

    Returns:
        dict: Formatted condition data.
    """
    compound_types = ["reagents", "catalyst", "solvent"]
    for compound_type in compound_types:
        name_key = f"{compound_type}_names"
        compound_name = conditions_data[name_key]
        for idx, name in enumerate(compound_name):
            if name == "Not Found":
                conditions_data[name_key][idx] = conditions_data[
                    f"{compound_type}_smiles"
                ][idx]
            elif name == "No Compound":
                conditions_data[name_key][idx] = "None"
    return conditions_data


def routes(
    solved_routes: Dict, selected_route: str, sustainability: Dict
) -> Tuple[dbc.Table, dbc.Table]:
    """
    Generates tables for route statistics and sustainability. For the routes table we only use the top condition set.

    Args:
        solved_routes: All the routes for a retrosynthesis
        selected_route: The route selected from the dropdown
        sustainability: The sustainability data which is used to colour the cell backgrounds

    Returns:
        tuple: Tables for route statistics and sustainability.
    """
    current_route = solved_routes[selected_route]
    # makes table with number of steps, route score, etc.
    route_sustainability = sustainability[selected_route]
    route_table = make_route_stats_table(
        current_route, selected_route, len(route_sustainability["steps"])
    )
    # make table headers based on number of steps, + 1 to use base 1 index
    table_headings = [
        html.Th(x + 1, style={"border-color": "black"})
        for x in range(len(route_sustainability["steps"]))
    ]
    table_headings.insert(0, html.Th("Step Analysis", style={"border-color": "black"}))
    table_header = [
        html.Thead(html.Tr(table_headings, style={"border-color": "black"}))
    ]
    # get table cell styles
    processed_sustainability = []
    for node_key, step_sustainability_dict in route_sustainability["steps"].items():
        conditions_dict = step_sustainability_dict.get("Condition Set 1")
        # update the sustainability dictionary to make the conditions table for the top condition set for each step.
        get_table_cell_styles(conditions_dict)
        conditions_dict.update({"node_id": node_key})
        processed_sustainability.append(conditions_dict)
    processed_sustainability.reverse()
    # make table rows
    row_keys = [
        "solvent",
        "temperature",
        "catalyst",
        "element_sustainability",
        "atom_economy",
        "safety",
        "weighted_median",
    ]
    table_rows = []
    for row_key in row_keys:
        title = row_key
        if title == "catalyst":
            title = "Stoichiometry/Catalyst"
        new_row = [
            html.Td(title.replace("_", " ").title(), style={"border-color": "black"})
        ]
        for step_dict in processed_sustainability:
            style = step_dict[row_key]["style"]
            cell = html.Td("", style=style)
            new_row.append(cell)
        table_rows.append(html.Tr(new_row, style={"border-color": "black"}))
    table_body = [html.Tbody(table_rows, style={"border-color": "black"})]
    sustainability_table = dbc.Table(
        table_header + table_body,
        bordered=True,
        style={"border-color": "black"},
        id="route-sustainability-table",
    )
    return route_table, sustainability_table


def make_route_stats_table(
    route: Dict, selected_route: str, number_of_steps: int
) -> dbc.Table:
    """
    Generates a table for route statistics.

    Args:
        route: ...
        selected_route: ...
        number_of_steps: ...

    Returns:
        dbc.Table: Table with route statistics.
    """
    route["score"] = round(route["score"], 2)
    df = pd.DataFrame(
        {
            "": ["Target Smiles", "Route", "Score", "Number of Steps"],
            "Route": [
                route["steps"][0]["smiles"],
                selected_route,
                route["score"],
                number_of_steps,
            ],
        }
    )
    route_table = dbc.Table.from_dataframe(
        df, bordered=True, hover=True, id="route-data-table"
    )
    return route_table


def generate_html_row(label: str, value: str, style: Optional[Dict] = None) -> html.Tr:
    """
    Generate an HTML row with a label and value.

    Args:
        label (str): The label for the row.
        value (Any): The value to be displayed in the row.
        style (Optional[Dict[str, Any]]): The optional styling for the value.

    Returns:
        dash_html_components.Tr: The HTML row component.
    """
    default_style = {
        "border-color": "black",
        "word-wrap": "break-word",
        "white-space": "normal",
    }
    if style is None:
        style = default_style
    return html.Tr(
        [
            html.Td(label, style=default_style),
            html.Td(value, style=style),
        ]
    )
