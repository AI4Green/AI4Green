from typing import Dict, List, Optional, Tuple, Union

from dash import html
from flask_login import current_user
from sources import services
from sources.auxiliary import get_workgroups


def make_workbooks_dropdown_options() -> Optional[Tuple[List[Dict], int]]:
    """
    Makes the workbooks dropdown with all workbooks the active user belongs to
    Returns:
        A list of dictionaries with workbook names as 'label' and primary key integer id as 'value'
        An integer of the primary key id of the first workbook in the list to act as a default value.

    """
    workgroups = get_workgroups()
    workbooks = [
        services.workbook.workbooks_from_workgroup(wg_name) for wg_name in workgroups
    ]
    # flatten nested lists
    workbooks = [workbook for sublist in workbooks for workbook in sublist]
    workbook_names = [wb.name for wb in workbooks]
    workbook_ids = [wb.id for wb in workbooks]
    workbook_options = [
        {"label": name, "value": int_id}
        for name, int_id in zip(workbook_names, workbook_ids)
    ]
    if workbook_ids:
        return workbook_options, workbook_ids[0]
    # if user is not in any workbooks
    return workbook_options, None


def make_conditions_dropdown(
    route_label: str,
    conditions_data: Dict[str, dict],
    weighted_sustainability_data: Dict[str, dict],
    tapped_node: Dict[str, Union[str, list]],
) -> Optional[Tuple[str, str, List[Dict]]]:
    """
    Generate dropdown options based on reaction conditions and sustainability.

    Args:
        route_label (str): The route information.
        conditions_data: Data containing reaction conditions.
        weighted_sustainability_data: Data containing sustainability information.
        tapped_node: The selected node information.

    Returns:
        Reaction conditions,
        Reaction sustainability,
        dropdown options.
        Returns (None, None, None) if condition prediction is unsuccessful.
    """
    # get data for the node the user has clicked.
    node_id = tapped_node["id"]
    conditions = conditions_data[route_label][node_id]
    sustainability = weighted_sustainability_data[route_label]["steps"][node_id]

    # return early if there is no condition data to show
    if conditions == "Condition Prediction Unsuccessful":
        return None, None, None

    # make the labels for the dropdown
    condition_set_labels = [f"Condition Set {x + 1}" for x in range(len(conditions))]

    # make the data and styling for the dropdown
    condition_set_styles = style_condition_set_dropdown(sustainability)

    condition_set_options = []
    for condition_set_label, style in zip(condition_set_labels, condition_set_styles):
        condition_set_options.append(
            {
                "label": html.Span([f"{condition_set_label}"], style=style),
                "value": condition_set_label,
                "style": {"width": "150%"},
                "width": "150%",
            }
        )
    return conditions, sustainability, condition_set_options


def style_condition_set_dropdown(reaction_weighted_sustainability: Dict) -> List[Dict]:
    """
    Creates the styling required to colour each option in the condition set dropdown by its sustainability
    Args:
        reaction_weighted_sustainability - dict with the weighted sustainability values for the reactions
    """
    hazard_colours = current_user.hazard_colors
    flag_score_to_phrase_dict = {
        4: "HighlyHazardous",
        3: "Hazardous",
        2: "Problematic",
        1: "Recommended",
    }
    # get weighted medians
    style_ls = []
    for condition_set in reaction_weighted_sustainability.values():
        weighted_median = round(condition_set["weighted_median"]["flag"], 0)
        chem21_phrase = flag_score_to_phrase_dict[weighted_median]
        background_colour = hazard_colours[chem21_phrase]
        text_colour = hazard_colours[f"{chem21_phrase}_text"]
        style_ls.append(
            {
                "background-color": background_colour,
                "color": text_colour,
                "width": "250%",
            }
        )
    return style_ls


def routes(
    active_retrosynthesis: Dict, active_weighted_sustainability: Dict
) -> List[Dict]:
    """
    Makes the routes dropdown where the user can change the active routes
    Args:
        active_retrosynthesis - the retrosynthesis currently being looked at contains all route labels
        active_weighted_sustainability - the weighted sustainability of the current retrosynthesis.
    Returns:
          The list of styled route dropdown options
    """
    number_of_steps_ls = [
        len(route["steps"])
        for route in active_weighted_sustainability["routes"].values()
    ]
    route_ls = [f"Route {x + 1}" for x in range(len(active_retrosynthesis["routes"]))]
    route_style_ls = style_routes_dropdown(active_weighted_sustainability)
    route_options = []
    for route_label, style, steps in zip(route_ls, route_style_ls, number_of_steps_ls):
        route_options.append(
            {
                "label": html.Span([f"{route_label} ({steps})"], style=style),
                "value": route_label,
                "style": {"width": "150%"},
                "width": "150%",
            }
        )
    return route_options


def style_routes_dropdown(weighted_sustainability: Dict) -> List[Dict]:
    """
    Style routes dropdown based on weighted sustainability scores.

    Args:
        weighted_sustainability Weighted sustainability data.

    Returns:
        A list of styles for routes dropdown based on weighted sustainability scores. 1 style per route.
    """
    hazard_colours = current_user.hazard_colors
    flag_score_to_phrase_dict = {
        4: "HighlyHazardous",
        3: "Hazardous",
        2: "Problematic",
        1: "Recommended",
    }
    # get weighted medians
    style_ls = []
    for route in weighted_sustainability["routes"].values():
        weighted_median = route["route_average"]["weighted_median"]
        weighted_median = round(weighted_median, 0)
        chem21_phrase = flag_score_to_phrase_dict[weighted_median]
        background_colour = hazard_colours[chem21_phrase]
        text_colour = hazard_colours[f"{chem21_phrase}_text"]
        style_ls.append(
            {
                "background-color": background_colour,
                "color": text_colour,
                "width": "250%",
            }
        )
    # convert weighted median number to hazard color
    return style_ls
