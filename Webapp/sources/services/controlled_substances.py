from sources import services, models
from sources.extensions import db
from typing import List, Optional
from flask import current_app
import os

def check_reaction_for_controlled_substances(reaction: models.Reaction):
    """
    Checks all chemicals in a reaction against AI4Green's list of controlled substances
    """
    substance_smiles = [
        reaction.reactants,
        reaction.reagents,
        reaction.products,
        reaction.solvent
    ]

    # this list of lists maintains the structure of reagents, reactants, products and solvents incase this is useful later
    checks = [check_smiles_in_controlled_substances(smi_list) for smi_list in substance_smiles]

    # if no compounds found return None
    if all(not sublist for sublist in checks):
        return None

    unique_structures = [
        substance
        for sublist in checks  # flatten the list of lists
        for substance in sublist
        if not check_duplicate(reaction, substance)  # filter out duplicates
    ]

    for substance in unique_structures:
        # add unique usages to db
        add(reaction, substance)

    return checks


def check_smiles_in_controlled_substances(smiles_list: Optional[List[str]]):
    query_inchis = [services.all_compounds.smiles_to_inchi(smi) for smi in smiles_list]
    return [inchi for inchi in query_inchis if inchi in current_app.config["CONTROLLED_SUBSTANCES"]]


def check_duplicate(reaction: models.Reaction, smiles: str) -> bool:
    duplicate = (
        db.session.query(models.ControlledSubstanceUsage)
        .filter(models.ControlledSubstanceUsage.reaction == reaction)
        .filter(models.ControlledSubstanceUsage.smiles == smiles)
    )

    if duplicate:
        return True

    return False


def add(reaction: models.Reaction, smiles: str) -> None:
    models.ControlledSubstanceUsage.create(
        creator=reaction.creator,
        workgroup=reaction.workgroup,
        workbook=reaction.workbook,
        reaction=reaction,
        controlled_substance_name=services.all_compounds.name_from_smiles(smiles),
        controlled_substance_smiles=smiles,
        controlled_substance_cas=services.all_compounds.cas_from_smiles(smiles),
    )


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

# Add controlled substance list to configs
current_app.config["CONTROLLED_SUBSTANCES"] = controlled_substance_inchi()
