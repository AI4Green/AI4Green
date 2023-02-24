from rdkit.Chem.Draw import rdMolDraw2D
from sources import db
from pony.orm import select, desc
from rdkit.Chem import AllChem
from rdkit import Chem
from operator import itemgetter
import re


def get_scheme_list(workbook, workgroup, sort_crit, size):
    # get reactions from reaction database
    if sort_crit == "time":
        reaction_list = list(select(reaction for reaction in db.Reaction if reaction.workbooks.name == workbook and
                                    reaction.workbooks.group.name == workgroup and reaction.status == "active").
                             order_by(desc(db.Reaction.time_of_creation))[:])
    else:
        reaction_list = list(select(reaction for reaction in db.Reaction if reaction.workbooks.name == workbook and
                                    reaction.workbooks.group.name == workgroup and reaction.status == "active").
                             order_by(db.Reaction.name)[:])
    scheme_list = []
    # get rxn schemes
    for j in range(len(reaction_list)):
        rxn_string = reaction_list[j].reaction_smiles
        # if the reaction string contains any letters - required for SMILES. Prevents '>>' breaking rdkit
        if re.search('[a-zA-Z]', rxn_string):
            # test to see if ions are present in which case further logic is needed
            rxn_string_split = rxn_string.split(' |')
            # if reaction scheme contains ions, longer method used in the if block is needed to make the scheme
            if len(rxn_string_split) > 1:
                # get list of ion positions
                ion_list = rxn_string.split("f:")[1].split('|')[0].split(',')
                # reagents to be added too.
                reactant_salt_dict_list, product_salt_dict_list = [], []
                # reactant_salt_list, reactant_salt_list_index, product_salt_list, product_salt_list_index = [], [], [], []
                # separate compounds, reactants, and products
                compounds_ls = rxn_string.split(' |')[0].replace('>>', '.').split('.')
                reactants_ls = rxn_string.split('>>')[0].split('.')
                products_ls = rxn_string.split('>>')[1].split(' |')[0].split('.')
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
                        insert_index = ion_indexes[0] - ions_in_previous_salts + previous_number_of_reactant_salts
                        reactant_salt_dict_list.append({"insert_index": insert_index, "end_index": insert_index + len(ion_indexes),
                                                       "number_of_ions": len(ion_indexes), "salt_smiles": salt})
                    else:
                        for salt_dict in product_salt_dict_list:
                            ions_in_previous_salts += salt_dict["number_of_ions"]
                        previous_number_of_product_salts = len(product_salt_dict_list)
                        insert_index = ion_indexes[0] - ions_in_previous_salts + previous_number_of_product_salts
                        product_salt_dict_list.append({"insert_index": insert_index, "end_index": insert_index + len(ion_indexes),
                                                       "number_of_ions": len(ion_indexes), "salt_smiles": salt})
                # replace the reactant ions with the new salt
                for salt_dict in reactant_salt_dict_list:
                    # eg [0, 1, 2, 3, 4] remove [0, 1] --> [2, 3, 4] --> insert 01 at index 0 --> [01, 2, 3, 4] -->
                    del reactants_ls[salt_dict["insert_index"]: salt_dict["end_index"]]
                    reactants_ls.insert(salt_dict["insert_index"], salt_dict["salt_smiles"])
                for salt_dict in product_salt_dict_list:
                    del products_ls[salt_dict["insert_index"]: salt_dict["end_index"]]
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
            scheme_list.append('')
    return scheme_list


def get_reaction_list(workbook, workgroup, sort_crit):
    # get reactions from reaction database
    reaction_list = list(select(reaction for reaction in db.Reaction if reaction.workbooks.name == workbook and
                                reaction.workbooks.group.name == workgroup and reaction.status == "active")[:])
    reactions = []
    for idx, reaction in enumerate(reaction_list):
        # for each reaction get the relevant info and shorten description if it's long
        description = reaction.description
        if reaction.creator.user:
            creator_email = reaction.creator.user.email
            creator_username = reaction.creator.user.username
        else:
            creator_email = "unknown"
            creator_username = "a deleted profile"
        if len(description) > 250:
            description = description[0:249] + "..."
        reaction_details = {'html_id': idx+1, 'name': reaction.name, 'description': description,
                            'time_of_creation': str(reaction.time_of_creation),
                            'reaction_smiles': reaction.reaction_smiles,
                            'reaction_table_data': reaction.reaction_table_data,
                            'summary_table_data': reaction.summary_table_data,
                            'workgroup': workgroup, 'workbook': workbook, 'completion_status': reaction.complete,
                            'reaction_id': reaction.reaction_id, 'creator_email': creator_email,
                            'creator_username': creator_username}
        reactions.append(reaction_details)
    if sort_crit == "time":
        reactions = sorted(reactions, key=itemgetter('time_of_creation'), reverse=True)
    else:
        reactions = sorted(reactions, key=itemgetter('name'))
    return reactions
