import re
from operator import itemgetter
from typing import Dict, List, Tuple

from flask import Response
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from sources import models
from sources.extensions import db


def get_scheme_list(
    workbook_name: str, workgroup_name: str, sort_crit: str, size: str
) -> List[str]:
    reaction_list = []
    query = (
        db.session.query(models.Reaction)
        .filter(models.Reaction.status == "active")
        .join(models.WorkBook)
        .filter(models.WorkBook.name == workbook_name)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup_name)
    )
    if sort_crit == "time":
        reaction_list = query.order_by(models.Reaction.time_of_creation.desc()).all()
    elif sort_crit == "AZ":
        reaction_list = query.order_by(models.Reaction.name.asc()).all()
    return make_scheme_list(reaction_list, size)


def make_scheme_list(reaction_list: List[models.Reaction], size: str) -> List[str]:
    scheme_list = []
    # get rxn schemes
    for j in range(len(reaction_list)):
        rxn_string = reaction_list[j].reaction_smiles
        # if the reaction string contains any letters - required for SMILES. Prevents '>>' breaking rdkit
        if re.search("[a-zA-Z]", rxn_string):
            # we test to see if ions are present in which case further logic is needed
            # first we see if it is from marvin js and contains ions
            if len(rxn_string.split(" |")) > 1:
                rxn = ion_containing_cx_smiles(rxn_string)
            elif "+" in rxn_string or "-" in rxn_string:
                rxn = ion_containing_smiles(rxn_string)
            # reactions with no ions - make rxn object directly from string
            else:
                rxn = AllChem.ReactionFromSmarts(rxn_string, useSmiles=True)
            # draw reaction scheme
            if size == "small":
                d2d = rdMolDraw2D.MolDraw2DSVG(400, 150)
            else:
                d2d = rdMolDraw2D.MolDraw2DSVG(600, 225)
            d2d.DrawReaction(rxn)
            # return drawing text
            scheme = d2d.GetDrawingText()
            scheme_list.append(scheme)
        else:
            scheme_list.append("")
    return scheme_list


def get_reaction_list(workbook, workgroup, sort_crit) -> List[Dict]:
    reaction_list = (
        db.session.query(models.Reaction)
        .filter(models.Reaction.status == "active")
        .join(models.WorkBook)
        .filter(models.WorkBook.name == workbook)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup)
        .all()
    )
    return process_reactions_to_dict(reaction_list, sort_crit)


def process_reactions_to_dict(reaction_list, sort_crit) -> List[Dict]:
    reactions = []
    for idx, reaction in enumerate(reaction_list):
        # for each reaction get the relevant info and shorten description if it's long
        description = reaction.description
        if reaction.creator_person.user:
            creator_email = reaction.creator_person.user.email
            creator_username = reaction.creator_person.user.username
        else:
            creator_email = "unknown"
            creator_username = "a deleted profile"
        if len(description) > 250:
            description = description[0:249] + "..."
        reaction_details = {
            "html_id": idx + 1,
            "name": reaction.name,
            "description": description,
            "time_of_creation": str(reaction.time_of_creation),
            "time_of_update": str(reaction.time_of_update),
            "reaction_smiles": reaction.reaction_smiles,
            "reaction_table_data": reaction.reaction_table_data,
            "summary_table_data": reaction.summary_table_data,
            "workgroup": reaction.WorkBook.WorkGroup.name,
            "workbook": reaction.WorkBook.name,
            "completion_status": reaction.complete,
            "reaction_id": reaction.reaction_id,
            "creator_email": creator_email,
            "creator_username": creator_username,
            "addenda": reaction.addenda,
        }
        reactions.append(reaction_details)
    return (
        sorted(reactions, key=itemgetter("time_of_creation"), reverse=True)
        if sort_crit == "time"
        else sorted(reactions, key=itemgetter("name"))
    )


def ion_containing_smiles(smiles: str) -> Chem.rdChemReactions.ChemicalReaction:
    """
    Process reaction smiles exported from ketcher that contain ions to a rdkit ChemicalReaction to make an image from
    """
    reactants, products = smiles.split(">>")
    reactants_list = reactants.split(".")
    ion_reactants = [
        (idx, reactant)
        for idx, reactant in enumerate(reactants_list)
        if "+" in reactants or "-" in reactants
    ]
    ionic_reactants = ions_to_ionic_compounds(ion_reactants)
    # from our index list we want to remove all of these from the reactants list and replace with a new ionic compound
    for idx, ionic_compound in enumerate(ionic_reactants):
        reactants_list.insert(
            ionic_compound["idx_list"][0] + idx, ionic_compound["string"]
        )
    # then remove each individual ion element from the list
    for ionic_compound in ionic_reactants:
        for ion in ionic_compound["component_ions"]:
            reactants_list.remove(ion)  # by elem identity not idx


def ions_to_ionic_compounds(ion_list: List[Tuple[int, str]]) -> List[Dict]:
    previous_ion_original_index = 99
    finished_ion_list = []
    current_ion = {"string": "", "idx_list": [], "component_ions": []}

    for idx, ion in enumerate(ion_list):
        # if ions have adjacent indexes, concatenate the strings
        if previous_ion_original_index == ion[0] - 1 or idx == 0:
            # adjacent ions should be combined
            current_ion["string"] += ion[1]
            current_ion["idx_list"].append(ion[0])
            current_ion["component_ions"].append(ion[1])
        # after 2 ions we check to see if this is done
        if len(current_ion["idx_list"]) >= 2:
            print("longer than 2", current_ion)
            # assess current ion for conditions of termination and adding current_ion_string to the finished list
            if assess_neutrality(current_ion["string"]) == "neutral":
                finished_ion_list.append(
                    current_ion["string"].replace("+", "").replace("-", "")
                )
                current_ion = {"string": "", "idx_list": [], "component_ions": []}
                continue
            # or if there are no adjacent ions we also reset. we can see this by seeing if the idx +1 is in a tup[le in ion_reactants
            if not [idx + 1 in ion_reactant_tuple for ion_reactant_tuple in ion_list]:
                # then we terminate again and assume the ionic compound is as is, as there are no more ions to add regardless of charge
                finished_ion_list.append(current_ion)
                current_ion = {"string": "", "idx_list": [], "component_ions": []}
                continue
        previous_ion_original_index = ion[0]
    return finished_ion_list


def assess_neutrality(ion_string: str) -> str:
    """
    We assess if the ion string is neutral and likely an ionic compound or if it is still charged and likely that more
    ions remain
    """
    # [Zn+2].[Cl-].[Cl-]>>O
    charge_balance = 0
    for idx, char in enumerate(ion_string):
        if char == "+" or char == "-":
            # if ion_string[idx + 1].isdigit():
            if len(ion_string) > idx + 1:
                charge_multiplier = (
                    ion_string[idx + 1] if ion_string[idx + 1].isdigit() else 1
                )
            else:
                charge_multiplier = 1
            if char == "+":
                charge_balance += charge_multiplier
            elif char == "-":
                charge_balance -= charge_multiplier
    return "neutral" if charge_balance == 0 else "not neutral"


def ion_containing_cx_smiles(cx_smiles: str) -> Chem.rdChemReactions.ChemicalReaction:
    """ "
    Process reaction smiles from marivn js that contain ions to a rdkit ChemicalReaction to make an image from
    These use CXSMILES to indicate any ions and their posiitons f:idx,idx2|SMILES
    """
    # get list of ion positions
    ion_list = cx_smiles.split("f:")[1].split("|")[0].split(",")
    # reagents to be added too.
    reactant_salt_dict_list, product_salt_dict_list = [], []
    # separate compounds, reactants, and products
    compounds_ls = cx_smiles.split(" |")[0].replace(">>", ".").split(".")
    reactants_ls = cx_smiles.split(">>")[0].split(".")
    products_ls = cx_smiles.split(">>")[1].split(" |")[0].split(".")
    # for each ionic compound pair the constituent ions to make a salt
    for idx, ion_group in enumerate(ion_list):
        ion_indexes = ion_group.split(".")
        ion_indexes = [int(x) for x in ion_indexes]
        salt = ".".join([compounds_ls[x] for x in ion_indexes])
        ions_in_previous_salts = 0
        if len(reactants_ls) > ion_indexes[-1]:
            # append salt to salt list, and compute+append the index of the ions which the salt will substitute
            for salt_dict in reactant_salt_dict_list:
                ions_in_previous_salts += salt_dict["number_of_ions"]
            previous_number_of_reactant_salts = len(reactant_salt_dict_list)
            insert_index = (
                ion_indexes[0]
                - ions_in_previous_salts
                + previous_number_of_reactant_salts
            )
            reactant_salt_dict_list.append(
                {
                    "insert_index": insert_index,
                    "end_index": insert_index + len(ion_indexes),
                    "number_of_ions": len(ion_indexes),
                    "salt_smiles": salt,
                }
            )
        else:
            for salt_dict in product_salt_dict_list:
                ions_in_previous_salts += salt_dict["number_of_ions"]
            previous_number_of_product_salts = len(product_salt_dict_list)
            insert_index = (
                ion_indexes[0]
                - ions_in_previous_salts
                + previous_number_of_product_salts
            )
            product_salt_dict_list.append(
                {
                    "insert_index": insert_index,
                    "end_index": insert_index + len(ion_indexes),
                    "number_of_ions": len(ion_indexes),
                    "salt_smiles": salt,
                }
            )
    # replace the reactant ions with the new salt
    for salt_dict in reactant_salt_dict_list:
        # eg [0, 1, 2, 3, 4] remove [0, 1] --> [2, 3, 4] --> insert 01 at index 0 --> [01, 2, 3, 4] -->
        del reactants_ls[salt_dict["insert_index"] : salt_dict["end_index"]]
        reactants_ls.insert(salt_dict["insert_index"], salt_dict["salt_smiles"])
    for salt_dict in product_salt_dict_list:
        del products_ls[salt_dict["insert_index"] : salt_dict["end_index"]]
        products_ls.insert(salt_dict["insert_index"], salt_dict["salt_smiles"])
    # convert to mol objects
    reactants_mols_ls = [Chem.MolFromSmiles(x) for x in reactants_ls]
    products_mols_ls = [Chem.MolFromSmiles(x) for x in products_ls]
    # reagents
    # make a rdkit Reaction object.
    rxn = Chem.rdChemReactions.ChemicalReaction()
    # add reactants, reagents, and products to the Reaction object
    for reactant in reactants_mols_ls:
        rxn.AddReactantTemplate(reactant)
    for product in products_mols_ls:
        rxn.AddProductTemplate(product)
    return rxn
