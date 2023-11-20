import re
from operator import itemgetter
from typing import Callable, Dict, List, Tuple

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

    Args:
        smiles (str): Reaction SMILES string.

    Returns:
        rdChemReactions.ChemicalReaction: rdkit ChemicalReaction object.
    """
    reactants, products = smiles.split(">>")
    reactants_list = update_ion_containing_list(reactants.split("."))
    products_list = update_ion_containing_list(products.split("."))

    reactants_mols_ls = [Chem.MolFromSmiles(x) for x in reactants_list]
    products_mols_ls = [Chem.MolFromSmiles(x) for x in products_list]

    if None in reactants_mols_ls or None in products_mols_ls:
        print("STOP")

    rxn = Chem.rdChemReactions.ChemicalReaction()
    add_mols_to_reaction(reactants_mols_ls, rxn.AddReactantTemplate)
    add_mols_to_reaction(products_mols_ls, rxn.AddProductTemplate)

    return rxn


def add_mols_to_reaction(mols: List[Chem.Mol], add_function: Callable) -> None:
    """
    Add molecules to a ChemicalReaction object.

    Args:
        mols (List[Union[Chem.Mol, None]]): List of rdkit Mol objects.
        add_function

    Returns:
        None
    """
    for mol in mols:
        if mol is not None:
            add_function(mol)


def update_ion_containing_list(compounds_list: List[str]) -> List[str]:
    """
    Update a list of compounds containing ions.

    Args:
        compounds_list (List[str]): List of compounds.

    Returns:
        List[str]: Updated list of compounds.
    """
    ion_reactants = [
        (idx, reactant)
        for idx, reactant in enumerate(compounds_list)
        if "+" in reactant or "-" in reactant
    ]
    ionic_reactants = ions_to_ionic_compounds(ion_reactants)

    for idx, ionic_compound in enumerate(ionic_reactants):
        compounds_list.insert(
            ionic_compound["idx_list"][0] + idx, ionic_compound["string"]
        )

    for ionic_compound in ionic_reactants:
        for ion in ionic_compound["component_ions"]:
            compounds_list.remove(ion)

    return compounds_list


def ions_to_ionic_compounds(ion_list: List[Tuple[int, str]]) -> List[Dict]:
    """
    Convert a list of ions to ionic compounds.

    Args:
        ion_list (List[Tuple[int, str]]): List of tuples containing index and ion strings.

    Returns:
        List[Dict]: List of dictionaries representing ionic compounds.
    """
    previous_ion_index = 99
    finished_ion_list = []
    current_ion = {"string": "", "idx_list": [], "component_ions": []}

    for idx, ion in enumerate(ion_list):
        if previous_ion_index == ion[0] - 1 or len(current_ion["idx_list"]) == 0:
            current_ion["string"] += ion[1]
            current_ion["idx_list"].append(ion[0])
            current_ion["component_ions"].append(ion[1])

        if len(current_ion["idx_list"]) >= 2:
            if assess_neutrality(current_ion["string"]) == "neutral":
                current_ion["string"] = rearrange_ionic_compound(
                    current_ion["component_ions"]
                )
                finished_ion_list.append(current_ion)
                current_ion = {"string": "", "idx_list": [], "component_ions": []}
                continue

            if not [idx + 1 in ion_reactant_tuple for ion_reactant_tuple in ion_list]:
                current_ion["string"] = rearrange_ionic_compound(
                    current_ion["component_ions"]
                )
                finished_ion_list.append(current_ion)
                current_ion = {"string": "", "idx_list": [], "component_ions": []}
                continue
        previous_ion_index = ion[0]

    return finished_ion_list


def rearrange_ionic_compound(component_ions: List[str]) -> str:
    """
    Rearrange the ionic compound.

    Args:
        component_ions (List[str]): List of component ions.

    Returns:
        str: Rearranged ionic compound.
    """
    result_list = sorted(component_ions, key=lambda x: abs(get_charge(x)))
    max_charge_element = max(result_list, key=lambda x: abs(get_charge(x)))
    result_list.remove(max_charge_element)
    result_list.insert(len(result_list) // 2, max_charge_element)
    rearranged_ion = (
        "".join(result_list).replace("+", "").replace("-", "").replace(r"\d+", "")
    )
    return re.sub(r"\d+", "", rearranged_ion)


def get_charge(element: str) -> int:
    """
    Get the charge of an element.

    Args:
    - element (str): Element string.

    Returns:
    - int: Charge of the element.
    """
    match = re.search(r"(\+|\-)(\d+)", element)
    return int(match.group(2)) if match else 0


def assess_neutrality(ion_string: str) -> str:
    """
    Assess the neutrality of an ion string.

    Args:
    - ion_string (str): Ion string.

    Returns:
    - str: 'neutral' if the ion string is neutral, 'not neutral' otherwise.
    """
    charge_balance = 0
    for idx, char in enumerate(ion_string):
        if char == "+" or char == "-":
            charge_multiplier = (
                int(ion_string[idx + 1])
                if len(ion_string) > idx + 1 and ion_string[idx + 1].isdigit()
                else 1
            )
            charge_balance += charge_multiplier if char == "+" else -charge_multiplier

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
