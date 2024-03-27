"""
This module receives a reaction from Marvin JS as a
GET request and renders the reaction table template
"""

import re
from typing import Dict, List, Tuple, Union
from urllib.parse import quote
from urllib.request import urlopen

from flask import current_app, jsonify, render_template, request
from flask_login import login_required
from rdkit import Chem
from rdkit.Chem import Descriptors
from sources import models, services
from sources.auxiliary import abort_if_user_not_in_workbook, smiles_symbols
from sources.dto import ReactionNoteSchema

from . import reaction_table_bp

if not current_app.config["DEBUG"]:
    print("Trying to import STOUT")
    try:
        from STOUT import translate_forward

        print("STOUT successfully imported")
    except Exception:
        print("Failed to import STOUT")
        pass
else:
    print("Application in debug mode. Not importing STOUT")


# Processing data from Marvin JS and creating reaction table
@reaction_table_bp.route("/_process", methods=["GET"])
def process():
    # must be logged in
    """This function receives reagents and product from browser, finds
    their IUPAC names and molar weights in PubChem, forms the lists of
    reagents and solvents. and renders the reaction table"""

    # get user workbook
    demo = request.args.get("demo")
    reaction = None
    workbook = None
    if demo != "demo":
        workgroup = request.args.get("workgroup")
        workbook_name = request.args.get("workbook")
        workbook = services.workbook.get_workbook_from_group_book_name_combination(
            workgroup, workbook_name
        )
        abort_if_user_not_in_workbook(workgroup, workbook_name, workbook)
        reaction_id = request.args.get("reaction_id")
        reaction = services.reaction.get_from_reaction_id_and_workbook_id(
            reaction_id, workbook.id
        )
    # get the SMILES string of reactants and products from the args and then replace the symbols
    reactants0 = request.args.get("reactants", 0, type=str)
    products0 = request.args.get("products", 0, type=str)
    if not reactants0 and products0:
        return jsonify({"error": "Missing data!"})

    reactants_smiles_list, products_smiles_list = get_reactants_and_products_list(
        reactants0, products0
    )

    reactant_data = {
        "molecular_weight_list": [],
        "name_list": [],
        "hazard_list": [],
        "density_list": [],
        "primary_key_list": [],
    }
    product_data = {
        "molecular_weight_list": [],
        "name_list": [],
        "hazard_list": [],
        "density_list": [],
        "primary_key_list": [],
        "table_numbers": [],
    }

    # Find reactants in database then add data to the dictionary
    for idx, reactant_smiles in enumerate(reactants_smiles_list, 1):
        novel_compound = False  # false but change later if true
        mol = Chem.MolFromSmiles(reactant_smiles)
        if mol is None:
            return jsonify({"error": f"Cannot process Reactant {idx} structure"})
        inchi = Chem.MolToInchi(mol)
        reactant = services.compound.get_compound_from_inchi(inchi)

        # if no match, then check the workbook collection of novel compounds
        if reactant is None:
            if demo == "demo":  # if in demo mode don't search novel compounds
                return jsonify({"reactionTable": "demo", "novelCompound": ""})

            reactant = (
                services.novel_compound.get_novel_compound_from_inchi_and_workbook(
                    inchi, workbook
                )
            )
            novel_compound = True

            # if no match is found we inform the user the compound is not in the database
            if reactant is None:
                reactant_name = iupac_convert(reactant_smiles)
                # generate molweight
                reactant_mol_wt = round(mol_weight_generate(reactant_smiles), 2)
                novel_reactant_html = render_template(
                    "_novel_compound.html",
                    component="Reactant",
                    name=reactant_name,
                    number=idx,
                    mw=reactant_mol_wt,
                    smiles=reactant_smiles,
                )
                return jsonify(
                    {"reactionTable": novel_reactant_html, "novelCompound": True}
                )

        # now we have the compound/novel_compound object, we can get all the data and add to reactant_data dict
        get_compound_data(reactant_data, reactant, novel_compound)
    number_of_reactants = len(reactant_data["name_list"])

    # Find products in database then add data to the dictionary
    for idx, product_smiles in enumerate(products_smiles_list, 1):
        novel_compound = False  # false but change later if true
        mol = Chem.MolFromSmiles(product_smiles)
        if mol is None:
            return jsonify({"error": f"Cannot process product {idx} structure"})
        inchi = Chem.MolToInchi(mol)
        product = services.compound.get_compound_from_inchi(inchi)

        # if no match, then check the workbook collection of novel compounds
        if product is None:
            if demo == "demo":  # if in demo mode don't search novel compounds
                return jsonify({"reactionTable": "demo", "novelCompound": ""})

            product = (
                services.novel_compound.get_novel_compound_from_inchi_and_workbook(
                    inchi, workbook
                )
            )
            novel_compound = True

            # if no match is found we inform the user the compound is not in the database
            if product is None:
                product_name = iupac_convert(product_smiles)
                # generate molweight
                product_mol_wt = round(mol_weight_generate(product_smiles), 2)
                novel_product_html = render_template(
                    "_novel_compound.html",
                    component="Product",
                    name=product_name,
                    number=idx,
                    mw=product_mol_wt,
                    smiles=product_smiles,
                )
                return jsonify(
                    {"reactionTable": novel_product_html, "novelCompound": True}
                )

        # now we have the compound/novel_compound object, we can get all the data and add to product_data dict
        get_compound_data(product_data, product, novel_compound)
        product_data["table_numbers"].append(number_of_reactants + idx)
    number_of_products = len(product_data["name_list"])

    # Reagents - There are too many for a unselected dropdown. Could do recently used within workbook
    # identifiers = reagent_name + reagent_cas
    identifiers = []

    # Solvents - keep solvents that are not novel compounds or are novel compounds within the current workbook
    if demo == "demo":
        sol_rows = services.solvent.get_default_list()
    else:
        sol_rows = services.solvent.get_workbook_list(workbook)

    r_class, r_classes = services.reaction_classification.classify_reaction(
        reactants_smiles_list, products_smiles_list
    )

    # Now it renders the reaction table template
    reaction_table = render_template(
        "_reaction_table.html",
        reactants=reactant_data["name_list"],
        reactant_mol_weights=reactant_data["molecular_weight_list"],
        reactant_densities=reactant_data["density_list"],
        reactant_hazards=reactant_data["hazard_list"],
        reactant_primary_keys=reactant_data["primary_key_list"],
        number_of_reactants=number_of_reactants,
        number_of_products=number_of_products,
        identifiers=identifiers,
        reactant_table_numbers=[],
        products=product_data["name_list"],
        product_mol_weights=product_data["molecular_weight_list"],
        product_hazards=product_data["hazard_list"],
        product_primary_keys=product_data["primary_key_list"],
        product_table_numbers=product_data["table_numbers"],
        reagent_table_numbers=[],
        reaction_table_data="",
        summary_table_data="",
        sol_rows=sol_rows,
        reaction=reaction,
        reaction_class=r_class,
        reaction_classes=r_classes,
    )
    return jsonify({"reactionTable": reaction_table})


def get_reactants_and_products_list(
    reactants: str, products: str
) -> Tuple[List[str], List[str]]:
    """
    Process reactants and products strings to obtain lists of reactant and product SMILES.
    Converts words into SMILES symbols and identifies the format as either CXSMILES or SMILES.
    If ions are present, it processes the ionic SMILES accordingly.

    Args:
        reactants (str): A string representing the reactants in the reaction.
        products (str): A string representing the products in the reaction.

    Returns:
        Tuple[List[str], List[str]]: A tuple containing lists of reactant and product SMILES.

    """
    reactants_smiles = smiles_symbols(reactants)
    products_smiles = smiles_symbols(products)

    # form the reaction_smiles. Just from reactands and products to exclude reagents/other data
    reaction_smiles = (reactants_smiles + ">>" + products_smiles).replace(",", ".")

    # we get the smiles straight from the sketcher to see if we have CXSmiles
    sketcher_smiles = smiles_symbols(request.args.get("reactionSmiles"))

    # [OH-].[Na+]>>[Cl-].[Cl-].[Zn++] |f:0.1,2.3.4|
    if "|" in sketcher_smiles:
        reaction_smiles += re.search(r" \|[^\|]*\|$", sketcher_smiles).group()
        (
            reactants_smiles_list,
            products_smiles_list,
        ) = services.ions.reactants_and_products_from_ionic_cx_smiles(reaction_smiles)
    elif "+" in reaction_smiles or "-" in reaction_smiles:
        (
            reactants_smiles_list,
            products_smiles_list,
        ) = services.ions.reactants_and_products_from_ionic_smiles(reaction_smiles)
        # reactions with no ions - make rxn object directly from string
    else:
        reactants_smiles_list, products_smiles_list = reactants_smiles.split(
            ","
        ), products_smiles.split(",")
    return reactants_smiles_list, products_smiles_list


def get_compound_data(
    compound_data: Dict,
    compound: Union[models.Compound, models.NovelCompound],
    novel_compound: bool,
):
    """
    Update compound data dictionary with information from the given compound object.

    Args:
        compound_data (Dict): A dictionary containing lists to store compound data.
        compound (Union[models.Compound, models.NovelCompound]): The compound or novel compound object.
        novel_compound (bool): A boolean flag indicating whether the compound is a novel compound.
    """

    # now we have the compound/novel_compound object, we can get all the data
    molecular_weight = (
        float(compound.molec_weight) if compound.molec_weight != "" else 0
    )
    compound_data["molecular_weight_list"].append(molecular_weight)

    compound_name = compound.name if compound.name != "" else "Not found"
    compound_data["name_list"].append(compound_name)

    compound_hazard = (
        compound.hphrase if compound.hphrase != "No hazard codes found" else "Unknown"
    )
    compound_data["hazard_list"].append(compound_hazard)

    compound_density = compound.density if compound.density != "" else "-"
    compound_data["density_list"].append(compound_density)

    if novel_compound:
        compound_data["primary_key_list"].append((compound.name, compound.workbook))
    else:
        compound_data["primary_key_list"].append(compound.id)


@reaction_table_bp.route("/_save_reaction_note", methods=["POST"])
@login_required
def save_reaction_note():
    """Saves an reaction_note to the reaction object"""
    workgroup_name = request.form["workgroup"]
    workbook_name = request.form["workbook"]
    workbook_object = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )
    reaction_id = request.form["reactionID"]
    reaction = services.reaction.get_from_reaction_id_and_workbook_id(
        reaction_id, workbook_object.id
    )
    reaction_note_text = request.form["reactionNoteText"]
    author = services.person.from_current_user_email()
    new_addendum = services.reaction.add_addendum(reaction, reaction_note_text, author)
    schema = ReactionNoteSchema()
    return jsonify({"reaction_note": schema.dump(new_addendum)})


def iupac_convert(smiles: str) -> str:
    """
    Tries to make the iupac name for a compound not in the database.
    First we try the CIR service. Second we try the STOUT-pypi python package
    """
    # print("Running CIR")
    try:
        url = (
            "http://cactus.nci.nih.gov/chemical/structure/"
            + quote(smiles)
            + "/iupac_name"
        )  # https://opsin.ch.cam.ac.uk/opsin/cyclopropane.png
        iupac_name = urlopen(url, [5]).read().decode("utf8")
        return iupac_name
    except Exception:
        print("failed CIR")
    print("trying STOUT")
    try:
        iupac_name = translate_forward(smiles)
        if iupac_name != "Could not generate IUPAC name for SMILES provided.":
            return iupac_name
    except Exception:
        print("STOUT failed")
    return ""


def mol_weight_generate(smiles: str) -> float:
    """
    Uses RDKit to calculate the molecular weight for a compound from its SMILES string

    Args:
        smiles - the SMILES of the compound of interest

    Returns:
        The molecular weight of the compound.
    """
    # MolWt accounts for the average across isotopes but ExactMolWt only takes the most abundant isotope.
    return Descriptors.MolWt(Chem.MolFromSmiles(smiles))
