from typing import Dict, List

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from sources import models, services
from sources.extensions import db


class ReactionSustainabilityFlags:
    """
    Class to get the sustainability flags for a reaction from the conditions
    """

    def __init__(self, conditions: Dict):
        """
        Creates a new instance of the ReactionSustainabilityFlags class
        Args:
            conditions - the conditions being assessed for their sustainability
        """
        self.conditions = conditions
        self.solvent_flag = 0
        self.temperature_flag = 0
        self.catalyst_flag = 0
        self.element_sustainability_flag = 0
        self.atom_economy_flag = 0
        self.safety_flag = 0
        self.unweighted_median = 0

    def to_dict(self) -> Dict[str, int]:
        """
        Creates a dictionary of the sustainability results
        Returns:
            Dictionary where the key is the metric value is the sustainability score.
        """
        return {
            "solvent": self.solvent_flag,
            "temperature": self.temperature_flag,
            "catalyst": self.catalyst_flag,
            "element_sustainability": self.element_sustainability_flag,
            "atom_economy": self.atom_economy_flag,
            "safety": self.safety_flag,
            "unweighted_median": self.unweighted_median,
        }

    def get_sustainability_flags(self) -> List:
        """
        Calls the relevant functions get the sustainability flag scores
        Returns - a list with each of the sustainability metrics
        """
        self.get_solvent_flag()
        self.get_temperature_flag()
        self.get_catalyst_flag()
        self.get_element_sustainability_flag()
        self.get_atom_economy_flag()
        self.get_safety_flag()
        self.get_unweighted_median()
        sustainability_flag_list = [
            self.solvent_flag,
            self.temperature_flag,
            self.catalyst_flag,
            self.element_sustainability_flag,
            self.atom_economy_flag,
            self.safety_flag,
        ]
        return sustainability_flag_list

    def get_unweighted_median(self):
        """Get the unweighted median for each sustainability metric"""
        self.unweighted_median = {
            "flag": np.median(
                [
                    self.solvent_flag["flag"],
                    self.temperature_flag["flag"],
                    self.catalyst_flag["flag"],
                    self.element_sustainability_flag["flag"],
                    self.atom_economy_flag["flag"],
                    self.safety_flag["flag"],
                ]
            )
        }

    def get_solvent_flag(self):
        """Assesses the solvent sustainability"""
        # Solvent flags - but reverse default order, so most hazardous is 4 and least is 1
        # flag rate dictionary. 5 is non chem21 and has average score of 2 - as is uncategorised.
        solvent_flag_rate = {1: 4, 2: 3, 3: 2, 4: 1, 5: 2}
        solvents = self.conditions["solvent_names"]
        solvent_flags = []  # solvent flag list
        for solvent_name in solvents:
            # no solvent used is sustainable and gives a score of 1.
            if solvent_name == "No solvent":
                solvent_flag_score = 1

            # if a solvent is used, try and find it in the database.
            else:
                solvent = services.solvent.chem21_solvent_from_name(solvent_name)
                # If a solvent has been found, take its flag score, if not assign non-chem21 value of 5.
                if isinstance(solvent, models.Solvent):
                    solvent_flag = solvent.flag
                else:
                    solvent_flag = 5
                solvent_flag_score = solvent_flag_rate[solvent_flag]

            solvent_flags.append(solvent_flag_score)
        self.solvent_flag = {
            "solvent": self.conditions["solvent_names"],
            "flag": np.mean(solvent_flags),
        }

    def get_temperature_flag(self):
        """Assesses the temperature sustainability"""
        conditions = self.conditions
        if 0 < conditions["temperature"] < 70:
            temp_flag = 1
        elif -20 < conditions["temperature"] < 140:
            temp_flag = 2
        else:
            temp_flag = 3
        self.temperature_flag = {
            "temperature": conditions["temperature"],
            "flag": temp_flag,
        }

    def get_catalyst_flag(self):
        """Assesses the catalyst sustainability"""
        conditions = self.conditions
        # green if catalyst is used or if no catalyst and no reagents
        if (
            conditions["catalyst"] == ""
            and conditions["reagents"] == ""
            or conditions["catalyst"] != ""
        ):
            catalyst_flag = 1
        else:
            # yellow if no catalyst used. Cannot assess stoichiometry of reagents, excess would be red
            catalyst_flag = 2
        self.catalyst_flag = {"catalyst": conditions["catalyst"], "flag": catalyst_flag}

    def get_element_sustainability_flag(self):
        """Assesses the element remaining supply sustainability"""
        # first get list of all smiles strings
        compound_keys = ["solvent", "reagents", "catalyst", "reactant", "product"]
        full_smiles_ls = []
        for compound_type in compound_keys:
            if not isinstance(self.conditions[compound_type], str):
                continue
            sub_smiles_ls = self.conditions[compound_type].split(".")
            if not isinstance(compound_type, str):
                continue
            if sub_smiles_ls != [""]:
                full_smiles_ls += sub_smiles_ls
        # get list of elements in reaction from the smiles list by converting to atoms and then chemical symbols
        element_symbols = set()
        for smiles in full_smiles_ls:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            for atom in mol.GetAtoms():
                symbol = atom.GetSymbol()
                element_symbols.add(symbol)
        # get sustainability colour and set element_sustainability variable according to least sustainable element present
        element_sustainability_set = set(
            y.colour
            for y in [
                db.session.query(models.Element.colour)
                .filter(models.Element.symbol == symbol)
                .first()
                for symbol in element_symbols
            ]
        )

        if "red" in element_sustainability_set:
            element_sustainability_flag = 3
            element_years = "5-50 years"
        elif "yellow" in element_sustainability_set:
            element_sustainability_flag = 2
            element_years = "50-500 years"
        else:
            element_sustainability_flag = 1
            element_years = "+500 years"
        self.element_sustainability_flag = {
            "years": element_years,
            "flag": element_sustainability_flag,
        }

    def get_atom_economy_flag(self):
        """Assesses the atom economy sustainability with some assumptions made"""
        # Requires equivalents assumption.
        # count product atoms - calculate mass.
        # count reactant atoms - calculate mass, assume stoichiometric reactant + reagent, and 5% cat. loading

        lhs_compound_keys = ["reactant", "reagents", "catalyst"]
        lhs_total_mass = 0
        for compound_type in lhs_compound_keys:
            if not isinstance(self.conditions[compound_type], str):
                continue
            sub_smiles_ls = self.conditions[compound_type].split(".")
            if sub_smiles_ls == [""]:
                continue
            for smiles in sub_smiles_ls:
                # smiles to molecular weight
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    continue
                molec_weight = Descriptors.MolWt(mol)
                if compound_type == "catalyst":
                    molec_weight = molec_weight * 0.05
                lhs_total_mass += molec_weight
        rhs_total_mass = 0
        for smiles in self.conditions["product"].split("."):
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            molec_weight = Descriptors.MolWt(mol)
            rhs_total_mass += molec_weight
        if rhs_total_mass == 0 or lhs_total_mass == 0:
            atom_economy = 100
        else:
            atom_economy = round(rhs_total_mass / lhs_total_mass * 100, 1)
        if atom_economy > 90:
            if atom_economy > 100:
                atom_economy = 100
            atom_economy_flag = 1
        elif atom_economy > 70:
            atom_economy_flag = 2
        else:
            atom_economy_flag = 3
        self.atom_economy_flag = {
            "atom_economy": atom_economy,
            "flag": atom_economy_flag,
        }

    def get_safety_flag(self):
        """Assesses the sustainability of the known hazard codes of the chemicals used in this reaction"""
        compounds = ["reactant", "product", "reagents", "solvent", "catalyst"]
        all_hazard_codes_list = []
        hazard_category_list = []
        for compound_type in compounds:
            compound_id_ls = self.conditions[f"{compound_type}_ids"]
            for compound_id in compound_id_ls:
                if compound_id in ["No Compound", "Not Found", None]:
                    continue
                compound_object = services.compound.get(compound_id)
                compound_hazard_codes = compound_object.hphrase
                if (
                    compound_hazard_codes is None
                    or compound_hazard_codes == "No hazard codes found"
                ):
                    continue
                compound_hazard_codes_list = compound_hazard_codes.split("-")
                all_hazard_codes_list += compound_hazard_codes_list
                # from here onwards this code may be altered to use chem21 categories
                for hazard_code in compound_hazard_codes_list:
                    # find category
                    hazard_code_record = (
                        db.session.query(models.HazardCode)
                        .filter(models.HazardCode.code == hazard_code)
                        .first()
                    )

                    hazard_category_list.append(hazard_code_record.category)

        if "VH" in hazard_category_list:
            safety_flag = 3
        elif "H" in hazard_category_list:
            safety_flag = 2
        else:
            safety_flag = 1
        all_hazard_codes_list = [
            x for x in all_hazard_codes_list if x != "No hazard codes found"
        ]
        self.safety_flag = {
            "hazard_code_list": all_hazard_codes_list,
            "flag": safety_flag,
        }


class RouteSustainabilityFlags:
    """
    Gets a route average using the TOP rated condition set for each reaction step.
    """

    def __init__(self, reactions: Dict[str, dict]):
        """
        Creates a new instance of the RouteSustainabilityFlags class
        Args:
            reactions - the reactions making up the route being assessed on its sustainability
        """
        self.reactions = reactions
        self.solvent_flags = []
        self.median_solvent_flag = 0
        self.temperature_flags = []
        self.median_temperature_flag = 0
        self.catalyst_flags = []
        self.median_catalyst_flag = 0
        self.element_sustainability_flags = []
        self.median_element_sustainability_flag = 0
        self.atom_economy_flags = []
        self.median_atom_economy_flag = 0
        self.safety_flags = []
        self.median_safety_flag = 0
        self.unweighted_median_sustainability_score = 0

    def get_sustainability_flags(self):
        """Gets the sustainability flags if conditions are present"""
        conditions_present = False
        for node_label, reaction_sustainability in self.reactions.items():
            if reaction_sustainability == "Condition Prediction Unsuccessful":
                continue
            conditions_present = True
            top_conditions = reaction_sustainability["Condition Set 1"]
            self.get_flags(top_conditions)
        if conditions_present is True:
            self.get_median_flags()

    def get_flags(self, top_conditions: Dict):
        """
        Gets the flag integer value for each metric
        Args:
            top_conditions - dictionary with the metrics for each of the top conditions for a specific transformation
        """
        self.solvent_flags.append(top_conditions["solvent"]["flag"])
        self.temperature_flags.append(top_conditions["temperature"]["flag"])
        self.catalyst_flags.append(top_conditions["catalyst"]["flag"])
        self.element_sustainability_flags.append(
            top_conditions["element_sustainability"]["flag"]
        )
        self.atom_economy_flags.append(top_conditions["atom_economy"]["flag"])
        self.safety_flags.append(top_conditions["safety"]["flag"])

    def get_median_flags(self):
        """Get median flags for each metric across the route and a route median flag score"""
        self.median_solvent_flag = np.median(self.solvent_flags)
        self.median_temperature_flag = np.median(self.temperature_flags)
        self.median_catalyst_flag = np.median(self.catalyst_flags)
        self.median_element_sustainability_flag = np.median(
            self.element_sustainability_flags
        )
        self.median_atom_economy_flag = np.median(self.atom_economy_flags)
        self.median_safety_flag = np.median(self.safety_flags)
        self.unweighted_median_sustainability_score = np.median(
            [
                self.median_solvent_flag,
                self.median_temperature_flag,
                self.median_catalyst_flag,
                self.median_element_sustainability_flag,
                self.median_atom_economy_flag,
                self.median_safety_flag,
            ]
        )

    def to_dict(self) -> dict:
        """
        Returns a dictionary with sustainability scores for each metric as a list for every reaction
        and a float as the median
        """
        solvent = {
            "solvent": {"list": self.solvent_flags, "median": self.median_solvent_flag}
        }
        temperature = {
            "temperature": {
                "list": self.temperature_flags,
                "median": self.median_temperature_flag,
            }
        }
        catalyst = {
            "catalyst": {
                "list": self.catalyst_flags,
                "median": self.median_catalyst_flag,
            }
        }
        element_sustainability = {
            "element_sustainability": {
                "list": self.element_sustainability_flags,
                "median": self.median_element_sustainability_flag,
            }
        }
        atom_economy = {
            "atom_economy": {
                "list": self.atom_economy_flags,
                "median": self.median_atom_economy_flag,
            }
        }
        safety = {
            "safety": {"list": self.safety_flags, "median": self.median_safety_flag}
        }
        dic = {
            **solvent,
            **temperature,
            **catalyst,
            **element_sustainability,
            **atom_economy,
            **safety,
        }
        return dic


class AllRouteSustainability:
    """Gets the sustainability scores for a series of routes"""

    def __init__(self, all_conditions: Dict[str, dict]):
        self.conditions = all_conditions
        self.route_sustainability = {}

    def get(self) -> Dict:
        """
        Makes a dictionary with the sustainability assessed for each step in a route
        Returns the route sustainability dictionary
        """
        for route_label, route in self.conditions.items():
            route_sustainability = {}
            # a node is a single reaction step
            for node_label, conditions in route.items():
                # skip failed examples
                if not conditions or conditions == "Condition Prediction Unsuccessful":
                    route_sustainability.update(
                        {node_label: "Condition Prediction Unsuccessful"}
                    )
                    continue

                reaction_step_sustainability = {}
                # assess sustainability for each condition set for a reaction and get the flag scores
                for condition_label, condition_set in conditions.items():
                    conditions_sustainability = self._get_condition_set_sustainability(
                        condition_set
                    )
                    reaction_step_sustainability.update(
                        {condition_label: conditions_sustainability}
                    )

                # update the dictionaries with the complete sustainability results
                route_sustainability.update({node_label: reaction_step_sustainability})

            route_flags = self._get_route_flags(route_sustainability)

            self.route_sustainability.update(
                {
                    route_label: {
                        "steps": route_sustainability,
                        "route_average": route_flags,
                    }
                }
            )
        return self.route_sustainability

    @staticmethod
    def _get_route_flags(route_sustainability: Dict[str, dict]) -> Dict[str, dict]:
        """Gets the route sustainability flags and returns as a dict"""
        route_sustainability_obj = RouteSustainabilityFlags(route_sustainability)
        route_sustainability_obj.get_sustainability_flags()
        return route_sustainability_obj.to_dict()

    @staticmethod
    def _get_condition_set_sustainability(
        condition_set: Dict[str, dict]
    ) -> Dict[str, int]:
        """Gets the condition set/reaction sustainability flags and returns as a dict"""
        sustainability = ReactionSustainabilityFlags(condition_set)
        sustainability.get_sustainability_flags()
        return sustainability.to_dict()


def weighted_median_for_route(
    route_sustainability: Dict, sustainability_weightings: List[int]
) -> float:
    """Gets a weighted median sustainability score for a route
    Args:
        route_sustainability - dictionary with sustainability data for the route
        sustainability_weightings - weighting for each metric based on list index. From frontend sliders

    Returns:
        median_median - float which is the weighted median of all the steps in the route
    """
    route_average = route_sustainability["route_average"]
    average_values = [x["median"] for x in route_average.values()]
    df = pd.DataFrame({"value": average_values, "weight": sustainability_weightings})
    median_median = get_weighted_median(df)
    return median_median


def weighted_median_for_each_step(
    route_sustainability: Dict[str, dict], sustainability_weightings: List[int]
):
    """
    route sustainability contains a dictionary for the sustainability of each step.
    Each condition set dictionary has the unweighted median replaced by the weighted median for that condition set
    """
    metric_list = [
        "solvent",
        "temperature",
        "catalyst",
        "element_sustainability",
        "atom_economy",
        "safety",
    ]
    for step_label, step_data in route_sustainability["steps"].items():
        if step_data == "Condition Prediction Unsuccessful":
            continue

        for condition_label, condition_set in step_data.items():
            # get the numerical flag for each metric
            metric_flag_list = [condition_set[metric]["flag"] for metric in metric_list]

            # from the dataframe get the weighted median of the flags and then update the original dictionary.
            df = pd.DataFrame(
                {"value": metric_flag_list, "weight": sustainability_weightings}
            )
            condition_set_weighted_median = get_weighted_median(df)
            del condition_set["unweighted_median"]
            condition_set["weighted_median"] = {"flag": None}
            condition_set["weighted_median"]["flag"] = condition_set_weighted_median


def get_weighted_median(df: pd.DataFrame) -> float:
    """
    Gets a weighted median for a step or a route
    Args:
        df - A dataframe with 2 columns. metric 'value' and 'weight'.

    Returns:
        The weighted median for the step or route

    """
    # sort
    df.sort_values("value", inplace=True)
    # find cumulative sum of weights
    cumsum = df["weight"].cumsum()
    # find 50th percentile weight
    cutoff = df["weight"].sum() / 2.0
    w_median = df["value"][cumsum >= cutoff].iloc[0]
    return w_median
