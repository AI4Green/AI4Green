import re
from operator import itemgetter
from typing import Any, Dict, List

from flask import Response
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

from sources import models
from sources.extensions import db


def get_scheme_list(workbook, workgroup, sort_crit, size, reaction_id) -> List[str]:
    reaction_list = []
    query = (
        db.session.query(models.Reaction)
        .filter(models.Reaction.status == "active")
        .join(models.WorkBook)
        .filter(models.WorkBook.name == workbook)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup)
    )
    if sort_crit == "time":
        reaction_list = query.order_by(models.Reaction.time_of_creation.desc()).all()
    elif sort_crit == 'AZ':
        reaction_list = query.order_by(models.Reaction.name.asc()).all()
    elif sort_crit == 'single':
        reaction_list = query.filter(models.Reaction.reaction_id == reaction_id).all()
    return make_scheme_list(reaction_list, size)


def make_scheme_list(reaction_list, size) -> List[str]:
    scheme_list = []
    # get rxn schemes
    for j in range(len(reaction_list)):
        rxn_string = reaction_list[j].reaction_smiles
        # if the reaction string contains any letters - required for SMILES. Prevents '>>' breaking rdkit
        if re.search("[a-zA-Z]", rxn_string):
            # test to see if ions are present in which case further logic is needed
            rxn_string_split = rxn_string.split(" |")
            # if reaction scheme contains ions, longer method used in the if block is needed to make the scheme
            if len(rxn_string_split) > 1:
                # get list of ion positions
                ion_list = rxn_string.split("f:")[1].split("|")[0].split(",")
                # reagents to be added too.
                reactant_salt_dict_list, product_salt_dict_list = [], []
                # reactant_salt_list, reactant_salt_list_index, product_salt_list, product_salt_list_index = [], [], [], []
                # separate compounds, reactants, and products
                compounds_ls = rxn_string.split(" |")[0].replace(">>", ".").split(".")
                reactants_ls = rxn_string.split(">>")[0].split(".")
                products_ls = rxn_string.split(">>")[1].split(" |")[0].split(".")
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
                    reactants_ls.insert(
                        salt_dict["insert_index"], salt_dict["salt_smiles"]
                    )
                for salt_dict in product_salt_dict_list:
                    del products_ls[salt_dict["insert_index"] : salt_dict["end_index"]]
                    products_ls.insert(
                        salt_dict["insert_index"], salt_dict["salt_smiles"]
                    )
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
                # for reagents
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
