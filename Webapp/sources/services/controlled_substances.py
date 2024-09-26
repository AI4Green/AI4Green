from sources import services, models
from typing import List, Optional
from flask import current_app
import os

def check_reaction_for_controlled_substances(reaction: models.Reaction):
    """
    Checks all chemicals in a reaction against AI4Green's list of controlled substances
    """
    check_reactants = check_smiles_in_controlled_substances(reaction.reactants)
    check_reagents = check_smiles_in_controlled_substances(reaction.reagents)
    check_products = check_smiles_in_controlled_substances(reaction.products)

    print(check_reactants, check_reagents, check_products)
    if not (check_reactants or check_reagents or check_products):
        return None

    else:
        pass



def check_smiles_in_controlled_substances(smiles_list: Optional[List[str]] = None):
    query_cas_numbers = [services.all_compounds.smiles_to_inchi(smi) for smi in smiles_list]
    return [x for x in query_cas_numbers if x in current_app.config["CONTROLLED_SUBSTANCES"]]


def controlled_substance_inchi():
    """Returns list of CAS numbers for controlled chemicals"""
    print("150 MISSING INCHI IN THE CONTROLLED SUBSTANCE LIST")
    with open(
            os.path.join(
                os.path.dirname(
                    os.path.dirname(
                        os.path.abspath(__file__)
                    )
                ), "static", "controlled_substances_inchi.txt"
            ), "r"
    ) as f:
        return set(f.read().splitlines())

current_app.config["CONTROLLED_SUBSTANCES"] = controlled_substance_inchi()
