"""
This module receives a reaction from Marvin JS as a
GET request and renders the reaction table template
"""

import re
from typing import Dict, List, Tuple, Union
from urllib.parse import quote
from urllib.request import urlopen

from flask import jsonify, render_template, request
from flask_login import login_required
from numpy import number
from rdkit import Chem
from sources import models, services
from sources.auxiliary import smiles_symbols
from sources.decorators import workbook_member_required
from sources.dto import ReactionNoteSchema

from . import reaction_table_bp


class SketcherCompound:
    """
    A compound that has been processed from the sketcher.
    This class should contain all the functions needed to process smiles from front end and include all the info needed to render the reaction table
    on a compound basis
    """

    def __init__(
        self,
        smiles,
        idx,
        polymer_indices,
        workbook,
        demo,
        number_of_reactants=0,
        is_product=False,
    ):
        self.smiles = smiles
        self.demo = demo
        self.workbook = workbook

        if is_product:
            self.reaction_component = "Product"
            self.idx = idx + number_of_reactants
        else:
            self.reaction_component = "Reactant"
            self.idx = idx

        self.is_polymer = False
        self.is_novel_compound = False
        self.novel_compound_table = None
        self.compound_data = {}

        self.errors = []

        self.check_polymer(polymer_indices)
        self.check_invalid_molecule()

        if not self.errors:
            self.process_compound()

    def process_compound(self):
        # find out if compound is a novel compound (if polymer then novel compound is always true)
        compound, novel_compound = get_compound_all_tables(
            self.smiles, self.workbook, self.is_polymer, self.demo
        )

        if compound is None:
            self.handle_new_novel_compound()

        else:
            self.is_novel_compound = novel_compound
            get_compound_data(self.compound_data, compound, novel_compound)

    def handle_new_novel_compound(self):
        if self.demo == "demo":
            self.errors.append(jsonify({"reactionTable": "Demo", "novelCompound": ""}))
            return

        compound_name = iupac_convert(self.smiles)
        # generate molweight
        mol_wt = services.all_compounds.mol_weight_from_smiles(self.smiles)
        novel_reactant_html = render_template(
            "_novel_compound.html",
            component=self.reaction_component,
            name=compound_name,
            number=self.idx,
            mw=mol_wt,
            smiles=self.smiles,
            polymer=self.is_polymer,
        )
        self.novel_compound_table = jsonify(
            {"reactionTable": novel_reactant_html, "novelCompound": True}
        )

    def check_polymer(self, polymer_indices):
        """
        Set is_polymer to True if idx in polymer indices and process smiles
        """
        if self.idx in polymer_indices:
            self.is_polymer = True
            self.smiles = services.polymer_novel_compound.find_canonical_repeat(
                self.smiles
            )
            self.check_polymer_dummy_atom()
            self.check_copolymer()

    def check_copolymer(self):
        if self.smiles.count("{+n}") > 1:
            self.errors.append(
                jsonify(
                    {
                        "error": f"Cannot process {self.reaction_component} {self.idx} structure: copolymers are not yet supported"
                    }
                )
            )

    def check_invalid_molecule(self):
        mol = Chem.MolFromSmiles(self.smiles)
        if mol is None:
            self.errors.append(
                jsonify(
                    {
                        "error": f"Cannot process {self.reaction_component} {self.idx} structure"
                    }
                )
            )

    def check_polymer_dummy_atom(self):
        if self.smiles == "":
            self.errors.append(
                jsonify(
                    {
                        "error": f"Cannot process Product {self.idx} structure: dummy atoms are not yet supported"
                    }
                )
            )


def check_compound_errors(compound_list: List[SketcherCompound]) -> Union[str, None]:
    """
    checks errors for lists of sketcher compounds
    """
    for compound in compound_list:
        print(compound.compound_data)
        if compound.errors:
            return compound.errors[0]
    return None


def check_novel_compounds(compound_list: List[SketcherCompound]) -> Union[str, None]:
    """
    should this fn be included in the errors fn?
    """
    for compound in compound_list:
        if compound.novel_compound_table:
            return compound.novel_compound_table
    return None


# Processing data from Marvin JS and creating reaction table
@reaction_table_bp.route("/_update_reaction_table", methods=["GET"])
@workbook_member_required
def update_reaction_table():
    # must be logged in
    """This function receives reagents and product from browser, finds
    their IUPAC names and molar weights in PubChem, forms the lists of
    reagents and solvents. and renders the reaction table

    We should refactor this function to provide each reagent with its associated data instead
    write clases for product and reactants, with like name

    """

    # get user workbook
    demo = request.args.get("demo")
    tutorial = request.args.get("tutorial")
    reaction = None
    workbook = None
    if demo != "demo" and tutorial != "yes":
        workgroup = request.args.get("workgroup")
        workbook_name = request.args.get("workbook")
        workbook = services.workbook.get_workbook_from_group_book_name_combination(
            workgroup, workbook_name
        )
        reaction_id = request.args.get("reaction_id")
        reaction = services.reaction.get_from_reaction_id_and_workbook_id(
            reaction_id, workbook.id
        )

    # get position of polymers in reaction or leave empty if none
    polymer_mode = request.args.get("polymer")
    polymer_indices = request.args.get("polymerIndices")
    if polymer_indices:
        polymer_indices = list(map(int, polymer_indices.split(",")))
    else:
        polymer_indices = list()

    # get the SMILES string of reactants and products from the args and then replace the symbols
    reactants0 = request.args.get("reactants", 0, type=str)
    products0 = request.args.get("products", 0, type=str)
    if not reactants0 and products0:
        return jsonify({"error": "Missing data!"})

    # keep
    reactants_smiles_list, products_smiles_list = get_reactants_and_products_list(
        reactants0, products0
    )

    reactants = [
        SketcherCompound(
            smiles=x,
            idx=idx + 1,
            polymer_indices=polymer_indices,
            workbook=workbook,
            demo=demo,
        )
        for idx, x in enumerate(reactants_smiles_list)
    ]
    number_of_reactants = len(reactants)

    products = [
        SketcherCompound(
            smiles=y,
            idx=idx + 1,
            polymer_indices=polymer_indices,
            workbook=workbook,
            demo=demo,
            number_of_reactants=number_of_reactants,
            is_product=True,
        )
        for idx, y in enumerate(products_smiles_list)
    ]
    number_of_products = len(products)

    # check errors first
    for compound_group in (reactants, products):
        # now handles co polymer, dummy atom and invalid molecule errors
        error = check_compound_errors(compound_group)
        if error:
            return error

        novel_compound_row = check_novel_compounds(compound_group)
        if novel_compound_row:
            return novel_compound_row

    # # Find reactants in database then add data to the dictionary
    # for idx, reactant_smiles in enumerate(reactants_smiles_list, 1):
    #     novel_compound = False  # false but change later if true
    #
    #     if idx in polymer_indices:
    #         polymer = True
    #         if reactant_smiles.count("{+n}") > 1:
    #             return jsonify(
    #                 {
    #                     "error": f"Cannot process Reactant {idx} structure: copolymers are not yet supported"
    #                 }
    #             )
    #         reactant_smiles = services.polymer_novel_compound.find_canonical_repeat(
    #             reactant_smiles
    #         )
    #
    #
    #     # check if valid SMILES
    #     mol = Chem.MolFromSmiles(reactant_smiles)
    #     if mol is None:
    #         return jsonify({"error": f"Cannot process Reactant {idx} structure"})
    #
    #     # find compound in database
    #     reactant, novel_compound = get_compound_all_tables(
    #         reactant_smiles, workbook, polymer, demo
    #     )
    #
    #     # if no match is found we inform the user the compound is not in the database
    #     if reactant is None:
    #         if demo == "demo":
    #             return jsonify({"reactionTable": "Demo", "novelCompound": ""})
    #
    #         reactant_name = iupac_convert(reactant_smiles)
    #         # generate molweight
    #         reactant_mol_wt = services.all_compounds.mol_weight_from_smiles(
    #             reactant_smiles
    #         )
    #         novel_reactant_html = render_template(
    #             "_novel_compound.html",
    #             component="Reactant",
    #             name=reactant_name,
    #             number=idx,
    #             mw=reactant_mol_wt,
    #             smiles=reactant_smiles,
    #             polymer=polymer,
    #         )
    #         return jsonify(
    #             {"reactionTable": novel_reactant_html, "novelCompound": True}
    #         )
    #
    #     # now we have the compound/novel_compound object, we can get all the data and add to reactant_data dict
    #     get_compound_data(reactant_data, reactant, novel_compound)
    # # keep
    # number_of_reactants = len(reactants)
    #
    # # Find products in database then add data to the dictionary
    # for idx, product_smiles in enumerate(products_smiles_list, 1):
    #     polymer = False
    #
    #     if (idx + number_of_reactants) in polymer_indices:  # if polymer:
    #         polymer = True
    #         if product_smiles.count("{+n}") > 1:
    #             return jsonify(
    #                 {
    #                     "error": f"Cannot process Product {idx} structure: copolymers are not yet supported"
    #                 }
    #             )
    #         product_smiles = services.polymer_novel_compound.find_canonical_repeat(
    #             product_smiles
    #         )
    #         if product_smiles == "":
    #             return jsonify(
    #                 {
    #                     "error": f"Cannot process Product {idx} structure: dummy atoms are not yet supported"
    #                 }
    #             )
    #
    #     # check if valid SMILES
    #     mol = Chem.MolFromSmiles(product_smiles)
    #     if mol is None:
    #         return jsonify({"error": f"Cannot process Product {idx} structure"})
    #
    #     # find compound in database
    #     product, novel_compound = get_compound_all_tables(
    #         product_smiles, workbook, polymer, demo
    #     )
    #
    #     # if no match is found we inform the user the compound is not in the database
    #     if product is None:
    #         if demo == "demo":
    #             return jsonify({"reactionTable": "Demo", "novelCompound": ""})
    #
    #         product_name = iupac_convert(product_smiles)
    #         # generate molweight
    #         product_mol_wt = services.all_compounds.mol_weight_from_smiles(
    #             product_smiles
    #         )
    #         novel_product_html = render_template(
    #             "_novel_compound.html",
    #             component="Product",
    #             name=product_name,
    #             number=idx,
    #             mw=product_mol_wt,
    #             smiles=product_smiles,
    #             polymer=polymer,
    #         )
    #         return jsonify({"reactionTable": novel_product_html, "novelCompound": True})
    #
    #     # now we have the compound/novel_compound object, we can get all the data and add to product_data dict
    #     get_compound_data(product_data, product, novel_compound)
    #     product_data["table_numbers"].append(number_of_reactants + idx)
    # number_of_products = len(product_data["name_list"])

    # Reagents - There are too many for a unselected dropdown. Could do recently used within workbook
    # identifiers = reagent_name + reagent_cas
    identifiers = []

    # Solvents - keep solvents that are not novel compounds or are novel compounds within the current workbook
    if demo == "demo":
        sol_rows = services.solvent.get_default_list()
    else:
        sol_rows = services.solvent.get_workbook_list(workbook)

    # change rxn table in polymer mode
    if polymer_mode.lower() == "true":
        reaction_table_html = "_polymer_reaction_table.html"
        r_class = None
        r_classes = None
    else:
        reaction_table_html = "_reaction_table.html"
        r_class, r_classes = services.reaction_classification.classify_reaction(
            reactants_smiles_list, products_smiles_list
        )

    # Now it renders the reaction table template
    reaction_table = render_template(
        reaction_table_html,
        reactants=reactants,
        number_of_products=number_of_products,
        identifiers=identifiers,
        reactant_table_numbers=[],
        products=products,
        # product_intended_dps=product_data["intended_dps"],
        reagent_table_numbers=[],
        reaction_table_data="",
        summary_table_data="",
        sol_rows=sol_rows,
        reaction=reaction,
        reaction_class=r_class,
        reaction_classes=r_classes,
        polymer_indices=polymer_indices,
    )
    return jsonify({"reactionTable": reaction_table})

    # # Now it renders the reaction table template
    # reaction_table = render_template(
    #     reaction_table_html,
    #     reactants=reactants,
    #     reactant_mol_weights=reactant_data["molecular_weight_list"],
    #     reactant_densities=reactant_data["density_list"],
    #     reactant_hazards=reactant_data["hazard_list"],
    #     reactant_primary_keys=reactant_data["primary_key_list"],
    #     number_of_reactants=number_of_reactants,
    #     number_of_products=number_of_products,
    #     identifiers=identifiers,
    #     reactant_table_numbers=[],
    #     products=products,
    #     product_mol_weights=product_data["molecular_weight_list"],
    #     product_hazards=product_data["hazard_list"],
    #     product_primary_keys=product_data["primary_key_list"],
    #     product_table_numbers=product_data["table_numbers"],
    #     # product_intended_dps=product_data["intended_dps"],
    #     reagent_table_numbers=[],
    #     reaction_table_data="",
    #     summary_table_data="",
    #     sol_rows=sol_rows,
    #     reaction=reaction,
    #     reaction_class=r_class,
    #     reaction_classes=r_classes,
    #     polymer_indices=polymer_indices,
    # )
    # return jsonify({"reactionTable": reaction_table})


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
    if "|f:" in sketcher_smiles:
        reaction_smiles += re.search(r" \|[^\|]*\|$", sketcher_smiles).group()
        (
            reactants_smiles_list,
            products_smiles_list,
        ) = services.ions.reactants_and_products_from_ionic_cx_smiles(reaction_smiles)
    elif re.findall(r"(?<!\{)\+", reaction_smiles) or re.findall(
        r"(?<!\{)-", reaction_smiles
    ):
        # find "+" or "-" in reaction_smiles but not {-} or {+n} for polymers
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


def get_compound_all_tables(smiles, workbook, polymer, demo):
    """
    Retrieves a compound from the database, checking the Compound, NovelCompound, and PolymerNovelCompound tables.

    Args:
        smiles (str): The SMILES string for the compound.
        workbook (str):
        polymer (bool): True if compound is a polymer.
        demo (str): "demo" if in demo mode.

    Returns:
        compound: The database object for the compound or None if it cannot be found.
        novel_compound (bool): True if compound not present in compound databases.


    """
    novel_compound = False  # false but change later if true
    if polymer:
        compound = None
    else:
        mol = Chem.MolFromSmiles(smiles)
        inchi = Chem.MolToInchi(mol)
        compound = services.compound.from_inchi(inchi)

    # if no match by inchi, then check the workbook collection of novel compounds
    if compound is None:
        if demo == "demo":  # if in demo mode don't search novel compounds
            return compound, True

        if polymer:
            compound = services.polymer_novel_compound.from_smiles_and_workbook(
                smiles, workbook
            )
        else:
            compound = services.novel_compound.from_inchi_and_workbook(inchi, workbook)
        novel_compound = True

    return compound, novel_compound


def get_compound_data(
    compound_data: Dict,
    compound: Union[models.Compound, models.NovelCompound, models.PolymerNovelCompound],
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
    compound_data["molecular_weight"] = molecular_weight

    compound_name = compound.name if compound.name != "" else "Not found"
    compound_data["name"] = compound_name

    compound_hazard = (
        compound.hphrase if compound.hphrase != "No hazard codes found" else "Unknown"
    )
    compound_data["hazards"] = compound_hazard

    compound_density = compound.density if compound.density != "" else "-"
    compound_data["density"] = compound_density

    if novel_compound:
        compound_data["primary_key"] = (compound.name, compound.workbook)
    else:
        compound_data["primary_key"] = compound.id


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
    First we try the CIR service.
    """
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
    return ""
