from typing import Dict, List, Literal, Optional, Union

import pandas as pd

from ..style_sheets import CytoscapeStyles as cytoStyle
from .utils import rdkit_smiles_to_image

node_border_dict = cytoStyle.node_border_dict


class RetrosynthesisCytoscape:
    """Class for making the cytoscape which is the interactive visual tool that displays the retrosynthesis"""

    def __init__(self, solved_routes: Dict, selected_route: str):
        """
        Creates a new instance of the RetrosynthesisCytoscape class
        Args:
            solved_routes - the dictionary containing the solved routes
            selected_route - the key used to access the selected route format e.g., 'Route 1'
        """

        self.solved_routes = solved_routes
        self.selected_route_idx = int(selected_route[-1]) - 1
        self.displayed_route = solved_routes[selected_route]
        self.node_df = self._make_node_dataframe()

    def make_cytoscape_elements(self):
        """Makes the nodes (molecules), and edges  which denote the relationships between molecules"""
        nodes = self.make_nodes()
        edges = self.make_edges()
        elements = nodes + edges
        return elements

    def make_cytoscape_stylesheet(self):
        """Makes the style sheet responsible for the molecular images, layout, and styling such as node borders"""
        chemical_image_stylesheet = self._make_chemical_image_stylesheet()
        new_stylesheet = cytoStyle.basic_stylesheet + chemical_image_stylesheet
        return new_stylesheet

    def _make_node_dataframe(self):
        """Makes a Pandas Dataframe with the node data from the dictionary"""
        # make dataframe where each row is a node, index is node_id, columns are smiles and node type
        node_ls = []
        for node in self.displayed_route["steps"]:
            node_type = self._get_node_type(node)
            node_reaction_smiles = self._get_node_reaction_smiles(node)
            node_ls.append(
                {
                    "smiles": node["smiles"],
                    "node_type": node_type,
                    "reaction_class": node["reaction_class"],
                    "reaction_smiles": node_reaction_smiles,
                }
            )
        node_ids = [x["node_id"] for x in self.displayed_route["steps"]]
        node_df = pd.DataFrame(node_ls, index=node_ids)
        return node_df

    @staticmethod
    def _get_node_reaction_smiles(node: Dict) -> Optional[str]:
        """
        Gets the reaction smiles to make a node if it has children (non-terminal)
        Args:
            node - the molecule node which may or may not be terminal
        Returns:
            the reaction smiles required to make the node if it is not terminal

        """
        reaction_smiles = ""
        if node["child_smiles"]:
            reactants = ".".join(node["child_smiles"])
            reaction_smiles = reactants + ">>" + node["smiles"]
        return reaction_smiles

    @staticmethod
    def _get_node_type(
        node: Dict,
    ) -> Union[Literal["terminal"], Literal["target"], Literal["normal"]]:
        """
        Gets the node type, whether it is terminal, entered target, or 'normal' if it has children but is not the target
        Args:
            node - the node we are getting the type of
        Returns:
            The node type. Either terminal, target, or normal
        """
        if not node["child_smiles"]:
            node_type = "terminal"
        elif node["node_id"] == "node-0":
            node_type = "target"
        else:
            node_type = "normal"
        return node_type

    def make_nodes(self) -> List[Dict]:
        """This makes the nodes with the corresponding id and smiles string as the label"""
        # make the coordinate lists
        x_ls, y_ls = [], []
        for node_id in self.node_df.index:
            depth = int(node_id.split("-")[1])
            row_position = int(node_id.split("-")[-1]) + 1
            number_of_elements_in_row = len(
                [x for x in self.node_df.index if int(x.split("-")[1]) == depth]
            )
            x = (
                500 / number_of_elements_in_row / 2 * row_position
            )  # if 1 elem 400/1 = 400 400/2=200
            y = depth * 50
            # print(f'{node_id=}, {x=}, {y=}')
            x_ls.append(x)
            y_ls.append(y)

        nodes = [
            {
                "data": {
                    "element": "node",
                    "id": node_id,
                    "smiles": smiles,
                    "reaction_smiles": reaction_smiles,
                    "label": reaction_class,
                    "type": node_type,
                },
                # 'position': {'x': x, 'y': y},
            }
            for node_id, smiles, reaction_smiles, reaction_class, node_type in zip(
                self.node_df.index,
                self.node_df["smiles"],
                self.node_df["reaction_smiles"],
                self.node_df["reaction_class"],
                self.node_df["node_type"],
            )
        ]
        return nodes

    def make_edges(self) -> List[Dict]:
        """This returns a list of edges to connect the nodes"""
        sources, targets, labels = [], [], []
        for step in self.displayed_route["steps"]:
            source = step["node_id"]
            # find children node id by their smiles
            for child in step["child_smiles"]:
                # find node with the smiles
                for idx, smiles in enumerate(self.node_df["smiles"]):
                    if child == smiles:
                        target_idx = idx
                target_id = self.displayed_route["steps"][target_idx]["node_id"]
                sources.append(source), targets.append(target_id), labels.append(
                    f'{step["smiles"]} to {child}'
                )
        edges = [
            {
                "data": {
                    "element": "edge",
                    "id": f"{source}_{target}",
                    "source": source,
                    "target": target,
                }
            }
            for source, target in zip(sources, targets)
        ]
        return edges

    def _make_chemical_image_stylesheet(self) -> List[Dict]:
        """This returns the stylesheet which matches the node id to the png string to show the molecule image"""
        img_ls = []
        for smiles in self.node_df["smiles"]:
            img_data = rdkit_smiles_to_image(smiles)
            img_ls.append(img_data)
        chemical_image_stylesheet = [
            {
                "selector": "#" + node_id,
                "style": {
                    "background-image": rxn_image,
                    "border-color": node_border_dict[node_type]["colour"],
                    "border-style": node_border_dict[node_type]["style"],
                    "border-width": node_border_dict[node_type]["width"],
                },
            }
            for node_id, rxn_image, node_type in zip(
                self.node_df.index,
                img_ls,
                self.node_df["node_type"],
            )
        ]
        return chemical_image_stylesheet
