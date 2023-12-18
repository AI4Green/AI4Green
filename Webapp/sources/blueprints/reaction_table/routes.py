"""
This module receives a reaction from Marvin JS as a
GET request and renders the reaction table template
"""
import ast
import json
from datetime import datetime
from urllib.parse import quote
from urllib.request import urlopen

from flask import abort, current_app, jsonify, render_template, request
from flask_login import current_user, login_required
from rdkit import Chem  # Used for converting smiles to inchi
from rdkit.Chem import Descriptors
from sources import models, services
from sources.auxiliary import abort_if_user_not_in_workbook, smiles_symbols
from sources.dto import ReactionNoteSchema

# render_template renders html templates
# request parses incoming request data and gives access to it
# jsonify is used to send a JSON response to the browser
from sources.extensions import db  # imports the module with auxiliary functions

from . import reaction_table_bp  # imports the blueprint of the reaction table route
from .reaction_classification import classify_reaction

if not current_app.config["DEBUG"]:
    try:
        from STOUT import translate_forward
    except Exception:
        pass


# Processing data from Marvin JS and creating reaction table
@reaction_table_bp.route("/_process", methods=["GET"])
def process():
    # must be logged in
    """This function receives reagents and product from browser, finds
    their IUPAC names and molar weights in PubChem, forms the lists of
    reagents and solvents. and renders the reaction table"""

    reaction_table_data = request.args.get("reaction_table")
    if reaction_table_data:
        reaction_table_data = ast.literal_eval(reaction_table_data)
    else:
        reaction_table_data = "no data"
    summary_table_data = request.args.get("summary_table")
    if summary_table_data and summary_table_data != "{}":
        summary_table_data = ast.literal_eval(summary_table_data)
    else:
        summary_table_data = "no data"
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
        reaction = (
            db.session.query(models.Reaction)
            .filter(models.Reaction.reaction_id == reaction_id)
            .join(models.WorkBook)
            .join(models.WorkGroup)
            .filter(models.WorkGroup.name == workgroup)
            .first()
        )

    # Reactants, reagents, and product
    reactants0 = request.args.get(
        "reactants", 0, type=str
    )  # gets the SMILES string of reactants
    products0 = request.args.get(
        "products", 0, type=str
    )  # gets the SMILES string of products
    reactant_mol_weights = []  # the list of reactant molar weights
    reactant_hazards = []  # the list of reactant hazard codes
    reactant_densities = []
    reactant_primary_keys = []
    reactant_table_numbers = []  # the list of reactant numbers in the reaction table
    reagent_table_numbers = []
    product_mol_weights = []  # the list of product molar weights
    product_hazards = []  # the list of product hazards
    product_primary_keys = []
    product_table_numbers = []  # the list of product numbers in the reaction table
    reactant_rdmols = []
    product_rdmols = []

    if reactants0 and products0:  # if reactants and products are received,
        # Reactants
        reactants = reactants0.split(
            ","
        )  # then the SMILES string is split into separate reactants
        i = 0  # index of the reactant list and their counter
        for reactant_smiles in reactants:  # goes through the reactant list
            # replaces 'minus/plus/sharp' with the desired symbol
            reactant_smiles = smiles_symbols(reactant_smiles)
            novel_compound = False  # false but change later if true
            # finds the first compound with a matching inchi
            mol = Chem.MolFromSmiles(reactant_smiles)  # convert smiles to mol
            reactant_rdmols.append(mol)
            try:
                reactant_inchi = Chem.MolToInchi(mol)  # convert mol to inchi
            except Exception:
                return jsonify({"error": f"Cannot process Reactant {i+1} structure"})
            reactant = (
                db.session.query(models.Compound)
                .filter(models.Compound.inchi == reactant_inchi)
                .first()
            )
            if (
                reactant is None
            ):  # if no match, then check the workbook collection of novel compounds
                if demo == "demo":  # if in demo mode don't search novel compounds
                    novel_reactant_html = "Demo"
                    return jsonify(
                        {"reactionTable": novel_reactant_html, "novelCompound": ""}
                    )
                reactant = (
                    db.session.query(models.NovelCompound)
                    .filter(models.NovelCompound.inchi == reactant_inchi)
                    .join(models.WorkBook)
                    .filter(models.WorkBook.id == workbook.id)
                    .first()
                )
                novel_compound = True
                if (
                    reactant is None
                ):  # if no match is found we inform the user the compound is not in the database
                    # generate name
                    reactant_name = iupac_convert(reactant_smiles)
                    # generate molweight
                    reactant_mol_wt = round(mol_weight_generate(reactant_smiles), 2)
                    novel_reactant_html = render_template(
                        "_novel_compound.html",
                        component="Reactant",
                        name=reactant_name,
                        number=i + 1,
                        mw=reactant_mol_wt,
                        smiles=reactant_smiles,
                    )
                    return jsonify(
                        {"reactionTable": novel_reactant_html, "novelCompound": True}
                    )
            reactant_mol_weight = reactant.molec_weight  # the reactant mol weight
            if reactant_mol_weight != "":  # checks if it is not empty
                reactant_mol_weight = float(reactant_mol_weight)  # converts it to float
            else:
                reactant_mol_weight = 0  # if it is empty, the zero value is assigned
            reactant_mol_weights.append(
                reactant_mol_weight
            )  # appends the reactant molar weight to their list

            reactant_name = reactant.name  # the reactant name
            if reactant_name == "":  # if the reactants name is empty,
                reactant_name = "Not found"  # then "Not found" is assigned
            reactants[i] = reactant_name  # otherwise, it's added to the reactant list

            reactant_hazard = reactant.hphrase
            if (
                reactant_hazard == "No hazard codes found"
            ):  # if the reactants name is empty,
                reactant_hazard = "Unknown"  # then "Not found" is assigned
            reactant_hazards.append(reactant_hazard)

            reactant_density = reactant.density
            if reactant_density == "":  # if the reactant density value is empty
                reactant_density = "-"
            reactant_densities.append(reactant_density)

            if novel_compound:
                # novel compound saved as name_book_composite to differentiate
                reactant_name = reactant.name
                reactant_wb = reactant.workbook
                reactant_pk = (reactant_name, reactant_wb)
            else:
                reactant_pk = reactant.id
            reactant_primary_keys.append(reactant_pk)

            i += 1  # increments the index and counter
        number_of_reactants = i

        # Products
        products = products0.split(
            ","
        )  # splits the SMILEs string into separate products
        k = 0  # index of the product list and their counter
        for product_smiles in products:  # goes through the product list
            novel_compound = False
            # replaces 'minus/plus/sharp' with the desired symbol
            product_smiles = smiles_symbols(product_smiles)
            mol = Chem.MolFromSmiles(product_smiles)  # convert smiles to mol
            product_rdmols.append(mol)
            try:
                product_inchi = Chem.MolToInchi(mol)  # convert mol to inchi
            except Exception:
                return jsonify({"error": f"Cannot process Product {k + 1} structure"})
            product = (
                db.session.query(models.Compound)
                .filter(models.Compound.inchi == product_inchi)
                .first()
            )
            if (
                product is None
            ):  # if no match, then check the workbook collection of novel compounds
                if demo == "demo":  # if in demo mode don't look
                    novel_product_html = "Demo"
                    return jsonify(
                        {"reactionTable": novel_product_html, "novelCompound": ""}
                    )
                product = (
                    db.session.query(models.NovelCompound)
                    .filter(models.NovelCompound.inchi == product_inchi)
                    .join(models.WorkBook)
                    .filter(models.WorkBook.id == workbook.id)
                    .first()
                )

                novel_compound = True
                if (
                    product is None
                ):  # if no match is found we inform the user the compound is not in the database
                    # if no match is found we inform the user the compound is not in the database
                    product_mol_wt = round(mol_weight_generate(product_smiles), 2)
                    product_name = iupac_convert(product_smiles)
                    novel_product_html = render_template(
                        "_novel_compound.html",
                        component="Product",
                        name=product_name,
                        number=k + 1,
                        mw=product_mol_wt,
                        smiles=product_smiles,
                    )
                    return jsonify(
                        {"reactionTable": novel_product_html, "novelCompound": ""}
                    )
            product_mol_weight = product.molec_weight  # the product mol weight
            if product_mol_weight != "":  # checks if it is not empty
                product_mol_weight = float(product_mol_weight)  # converts it to float
            else:
                product_mol_weight = 0  # if it is empty, the zero value is assigned
            product_mol_weights.append(
                product_mol_weight
            )  # appends the product molar weight to their list

            product_name = product.name  # the product name
            if product_name == "":  # if the product name is empty,
                product_name = "Not found"  # then "Not found" is assigned
            products[k] = product_name  # assigning product name for the reaction table

            product_hazard = product.hphrase
            if (
                product_hazard == "No hazard codes found"
            ):  # if the reactants name is empty,
                product_hazard = "Unknown"  # then "Not found" is assigned
            product_hazards.append(product_hazard)

            if novel_compound:
                # novel compound saved as name_book_composite to differentiate
                product_name = product.name
                product_wb = product.workbook
                product_pk = (product_name, product_wb)
            else:
                product_pk = product.id
            product_primary_keys.append(product_pk)

            k += 1  # increments the index and counter
            product_table_numbers.append(
                number_of_reactants + k
            )  # assigning the product number in the reaction table
        number_of_products = k

        # Reagents - There are too many for a unselected dropdown.
        # identifiers = reagent_name + reagent_cas
        identifiers = []
        # Solvents - keep solvents that are not novel compounds or are novel compounds within the current workbook
        sol_rows = db.session.query(models.Solvent).all()
        if demo == "demo":
            sol_rows = [x for x in sol_rows if x.novel_compound == []]
        else:
            sol_rows = [
                x
                for x in sol_rows
                if x.novel_compound == []
                or (
                    workbook is not None and x.novel_compound[0].workbook == workbook.id
                )
            ]
        r_class, r_classes = classify_reaction(reactant_rdmols, product_rdmols)
        # Now it renders the reaction table template
        reaction_table = render_template(
            "_reaction_table.html",
            reactants=reactants,
            reactant_mol_weights=reactant_mol_weights,
            reactant_densities=reactant_densities,
            reactant_hazards=reactant_hazards,
            reactant_primary_keys=reactant_primary_keys,
            number_of_reactants=number_of_reactants,
            number_of_products=number_of_products,
            identifiers=identifiers,
            reactant_table_numbers=reactant_table_numbers,
            products=products,
            product_mol_weights=product_mol_weights,
            product_hazards=product_hazards,
            product_primary_keys=product_primary_keys,
            product_table_numbers=product_table_numbers,
            reagent_table_numbers=reagent_table_numbers,
            reaction_table_data=json.dumps(reaction_table_data),
            summary_table_data=json.dumps(summary_table_data),
            sol_rows=sol_rows,
            reaction=reaction,
            reaction_class=r_class,
            reaction_classes=r_classes,
        )
        return jsonify({"reactionTable": reaction_table})

    return jsonify({"error": "Missing data!"})


@reaction_table_bp.route("/_save_reaction_note", methods=["POST"])
@login_required
def save_reaction_note():
    """Saves an reaction_note to the reaction object"""
    workgroup = request.form["workgroup"]
    workbook = request.form["workbook"]
    workbook_object = (
        db.session.query(models.WorkBook)
        .filter(models.WorkBook.name == workbook)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup)
        .first()
    )
    reaction_id = request.form["reactionID"]
    reaction = (
        db.session.query(models.Reaction)
        .filter(models.Reaction.reaction_id == reaction_id)
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook_object.id)
        .first()
    )

    reaction_note_text = request.form["reactionNoteText"]
    author = (
        db.session.query(models.Person)
        .join(models.User)
        .filter(models.User.email == current_user.email)
        .first()
    )
    new_addendum = models.ReactionNote(
        text=reaction_note_text,
        time_of_creation=datetime.now(),
        author=author.id,
        reaction=reaction.id,
    )
    db.session.add(new_addendum)
    db.session.commit()

    schema = ReactionNoteSchema()
    return jsonify({"reaction_note": schema.dump(new_addendum)})


def iupac_convert(ids):
    print("Running CIR")
    try:
        url = (
            "http://cactus.nci.nih.gov/chemical/structure/" + quote(ids) + "/iupac_name"
        )  # https://opsin.ch.cam.ac.uk/opsin/cyclopropane.png
        ans = urlopen(url, [5]).read().decode("utf8")
        return ans
    except Exception:
        print("failed CIR")
    print("trying STOUT")
    try:
        ans = translate_forward(ids)
        return ans
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
