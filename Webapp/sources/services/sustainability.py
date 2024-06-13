from typing import Dict, List, Optional, Set

import rdkit
from flask import request
from sources import db, models


class SustainabilityMetrics:
    def __init__(
        self, reactant_data: Dict, reagent_data: Dict, solvent_data, product_data: Dict
    ):
        self.reactant_data = reactant_data
        self.reagent_data = reagent_data
        self.solvent_data = solvent_data
        self.product_data = product_data
        self.sustainability_data = {}

    def get_sustainability_metrics(self) -> Dict:
        self.sustainability_data["ae"] = self.get_atom_economy()
        self.sustainability_data["ae_flag"] = metric_flag(
            self.sustainability_data["ae"]
        )
        self.sustainability_data[
            "element_sustainability"
        ] = self.get_element_sustainability()
        self.sustainability_data[
            "element_sustainability_flag"
        ] = self.get_element_sustainability_flag()
        self.sustainability_data["solvent_flags"] = self.get_solvent_flags()
        return self.sustainability_data

    def get_solvent_flags(self) -> List[int]:
        # Solvent flags
        flag_rate = {
            1: "hazard-highly-hazardous",
            2: "hazard-hazardous",
            3: "hazard-warning",
            4: "hazard-acceptable",
            5: "non-chem21",
        }  # flag rate dictionary
        solvent_flags = []  # solvent flag list

        if self.solvent_data["solvents"][0] and self.solvent_data.get("solvents") != [
            " "
        ]:
            for solvent in self.solvent_data["solvents"]:
                solvent_flag = (
                    db.session.query(models.Solvent.flag)
                    .filter(models.Solvent.name == solvent)
                    .first()
                )
                solvent_flag = solvent_flag[0] if solvent_flag else None
                try:
                    solvent_flags.append(
                        flag_rate[solvent_flag]
                    )  # appends solvent flag to their list
                except KeyError:
                    solvent_flags.append(5)
        return solvent_flags

    def get_element_sustainability(self, reaction_smiles: str = None):
        # reaction smiles initially only has product and reactant SMILES strings as r1.rn>>p1.pn
        if reaction_smiles is None:
            reaction_smiles = str(request.form["reactionSmiles"])
        full_reaction_smiles_list = self.make_reaction_smiles_list(reaction_smiles)
        element_symbols_set = get_element_set(full_reaction_smiles_list)
        element_sustainability_set = get_element_sustainability_set(element_symbols_set)
        return element_sustainability_from_set(element_sustainability_set)

    def get_element_sustainability_flag(self):
        element_flag_dict = {
            "5-50 years": "hazard-hazardous",
            "50-500 years": "hazard-warning",
            "+500 years": "hazard-acceptable",
            "-select-": "hazard-reset-hazard",
        }
        return element_flag_dict[self.sustainability_data["element_sustainability"]]

    def get_atom_economy(self) -> Optional[float]:
        return (
            round(
                100
                * float(
                    self.product_data["product_molecular_weights"][
                        self.product_data["main_product_index"]
                    ]
                )
                / (
                    self.reactant_data["reactant_molecular_weight_sum"]
                    + self.reagent_data["reagent_molecular_weight_sum"]
                ),
                1,
            )
            if (
                self.reactant_data["reactant_molecular_weight_sum"]
                + self.reagent_data["reagent_molecular_weight_sum"]
            )
            > 0
            else 0
        )

    def make_reaction_smiles_list(self, reaction_smiles: str) -> List:
        reaction_smiles_ls = reaction_smiles.replace(">>", ".").split(".")
        full_reaction_smiles_ls = [
            x
            for x in reaction_smiles_ls
            + self.reagent_data["reagent_smiles_ls"]
            + self.solvent_data["solvent_smiles_ls"]
            if x
        ]
        return [x for x in full_reaction_smiles_ls if x]


def element_sustainability_from_set(element_sustainability_set: Set[str]) -> str:
    if "red" in element_sustainability_set:
        element_sustainability = "5-50 years"
    elif "yellow" in element_sustainability_set:
        element_sustainability = "50-500 years"
    elif "lime" in element_sustainability_set:
        element_sustainability = "+500 years"
    else:
        element_sustainability = "-select-"
    return element_sustainability


def get_element_set(reaction_smiles_list: List):
    element_symbols = set()
    for component in reaction_smiles_list:
        mol = rdkit.Chem.MolFromSmiles(component)
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            element_symbols.add(symbol)
    return element_symbols


def get_element_sustainability_set(element_symbols: Set[str]) -> Set[str]:
    return set(
        y.colour
        for y in [
            db.session.query(models.Element.colour)
            .filter(models.Element.symbol == symbol)
            .first()
            for symbol in element_symbols
        ]
    )


def metric_flag(metric: int) -> str:
    """
    Colour flags for metrics
    """
    if metric > 90:
        return "hazard-acceptable"
    elif metric < 70:
        return "hazard-hazardous"
    else:
        return "hazard-warning"


def percent_yield_from_rxn_data(rxn_data: Dict, obtained_yield_mass: float) -> float:
    """
    Calculates the percentage yield from rxn_data from the attribute: models.reaction.reaction_data.
    Args:
        rxn_data: dictionary with limiting reactant and desired product data
        obtained_yield_mass: the product mass yield - obtained from models.reaction.summary_data['real_product_mass']

    Returns:
        The percentage yield relative to the theoretical yield for a reaction.
    """
    # unit dicts
    mass_factor_dict = {"μg": 0.000001, "mg": 0.001, "g": 1}
    amount_factor_dict = {"μmol": 0.000001, "mmol": 0.001, "mol": 1}
    # get the amount of limiting reactant standardised for different amount units
    limiting_reactant_index = int(rxn_data["limiting_reactant_table_number"]) - 1
    amount_of_limiting_reactant = float(
        rxn_data["reactant_amounts_raw"][limiting_reactant_index]
    )
    amount_factor = amount_factor_dict[rxn_data["amount_units"]]
    normalised_limiting_reactant_amount = amount_of_limiting_reactant * amount_factor
    # get the theoretical yield of product standardised for different mass units
    main_product_index = int(rxn_data["main_product"]) - 1
    desired_product_molecular_weight = float(
        rxn_data["product_molecular_weights"][main_product_index]
    )
    mass_factor = mass_factor_dict[rxn_data["product_mass_units"]]
    theoretical_yield = (
        normalised_limiting_reactant_amount * desired_product_molecular_weight
    )
    normalised_yield_mass = obtained_yield_mass * mass_factor
    return round(normalised_yield_mass / theoretical_yield * 100, 2)
