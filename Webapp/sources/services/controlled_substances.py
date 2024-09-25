from sources import services, models
from typing import List, Optional, Union
import os

def check_reaction_for_controlled_substances(reaction: models.Reaction):
    """
    Checks all chemicals in a reaction against AI4Green's list of controlled substances
    """
    check_reactants = check_smiles_in_controlled_substances(reaction.reactants)
    check_reagents = check_smiles_in_controlled_substances(reaction.reagents)
    check_products = check_smiles_in_controlled_substances(reaction.products)

    if not (check_reactants or check_reagents or check_products):
        return None

    else:
        pass



def check_smiles_in_controlled_substances(smiles_list: Optional[List[str]] = None):
    query_cas_numbers = [services.all_compounds.cas_from_smiles(smi) for smi in smiles_list]
    return [x for x in query_cas_numbers if x in controlled_substances()]


def controlled_substances():
    """Returns list of CAS numbers for controlled chemicals"""
    with open(
            os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "", "controlled_substances.txt"
        ), "r"
    ) as f:
        return f.read().splitlines()
