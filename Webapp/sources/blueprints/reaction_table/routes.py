"""
This module receives a reaction from Marvin JS as a
GET request and renders the reaction table template
"""

import re
from typing import Dict, List, Tuple, Union
from urllib.parse import quote
from urllib.request import urlopen

from flask import json, jsonify, render_template, request
from flask_login import login_required
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
        workbook,
        demo,
        reaction_component,
        reaction_component_idx,
        polymer_indices=None,
        reaction_smiles="",
        reload=False,
    ):
        self.smiles = smiles
        self.inchi = ""
        self.idx = idx
        self.reaction_component_idx = reaction_component_idx
        self.demo = demo
        self.workbook = workbook
        self.reaction_component = reaction_component
        self.is_novel_compound = False
        self.is_polymer = False
        self.novel_compound_table = None
        self.compound_data = {}
        self.reload = reload
        self.errors = []

        self.check_for_polymer(polymer_indices, reaction_smiles)

        self.check_invalid_molecule()
        self.check_polymer_dummy_atom()
        self.check_copolymer()

        if not self.errors and self.reaction_component != "Solvent":
            self.process_compound()

    def process_compound(self):
        # find out if compound is a novel compound (if polymer then novel compound is always true)# fails for polymers lmao
        compound, novel_compound = get_compound_all_tables(
            self.smiles, self.workbook, self.is_polymer, self.demo
        )
        # only check for novel compound if reaction is not being reloaded
        if compound is None:
            self.handle_new_novel_compound()

        else:
            self.is_novel_compound = novel_compound
            get_compound_data(self.compound_data, compound, novel_compound)

    def check_for_polymer(self, polymer_indices, reaction_smiles):
        if polymer_indices is not None:
            self.check_polymer_indices_for_polymer(polymer_indices)
        else:
            self.check_reaction_smiles_for_polymer(reaction_smiles)

    def check_polymer_indices_for_polymer(self, polymer_indices):
        """
        Set is_polymer to True if idx in polymer indices and process smiles
        """
        if self.idx in polymer_indices:
            self.is_polymer = True
            self.smiles = services.polymer_novel_compound.find_canonical_repeat(
                self.smiles
            )

    def check_reaction_smiles_for_polymer(self, reaction_smiles):
        reactant_smiles, product_smiles = get_reactants_and_products_list(
            reaction_smiles
        )
        smiles_list = (
            reactant_smiles if self.reaction_component == "Reactant" else product_smiles
        )
        # only check reactant and products for polymers
        if self.reaction_component in ["Reactant", "Product"]:
            if "{+n}" in smiles_list[self.reaction_component_idx]:
                self.is_polymer = True

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
            # chenage for novel compound table
            number=self.idx,
            mw=mol_wt,
            smiles=self.smiles,
            polymer=self.is_polymer,
        )
        self.novel_compound_table = jsonify(
            {"reactionTable": novel_reactant_html, "novelCompound": True}
        )

    def add_solvent_sustainability_flags(self):
        flag = services.solvent.sustainability_from_primary_key(
            self.compound_data["ids"]
        )
        self.compound_data[
            "sustainability_flag"
        ] = services.solvent.convert_sustainability_flag_to_text(flag)

    @classmethod
    def from_reaction_table_dict(cls, reaction_table_dict, workbook):
        """
        Create SketcherCompound instances from a dictionary of reaction data.

        Args:
            reaction_table_dict (dict): Dictionary containing reaction data with keys for reactants, products, and other details.
            workbook:

        Returns:
            dict: A list of SketcherCompound instances representing individual reactants and products.
        """

        component_lists = {"reactant": [], "reagent": [], "solvent": [], "product": []}
        units = {}

        reaction_param_keys = [
            "limiting_reactant_table_number",
            "main_product",
            "mass_units",
            "polymerisation_type",
            "amount_units",
            "volume_units",
            "solvent_volume_units",
            "product_mass_units",
            "product_amount_units",
        ]

        for param in reaction_param_keys:
            units[param] = reaction_table_dict.pop(param)

        number_of_compounds = 0
        for component_type, component_list in component_lists.items():
            # get all relevant items from reaction_table_dict
            sub_dict = {
                k: v for k, v in reaction_table_dict.items() if component_type in k
            }
            for idx, name in enumerate(sub_dict.get(component_type + "_names")):
                number_of_compounds += 1
                compound_data = {}
                for key in sub_dict.keys():
                    # new_key = key.replace(component_type + "_", "")
                    value = sub_dict.get(key, "")
                    if "units" in key:
                        compound_data[key.replace(component_type + "_", "")] = value
                    else:
                        try:
                            compound_data[
                                key.replace(component_type + "_", "")
                            ] = value[idx]
                        except IndexError:
                            compound_data[key.replace(component_type + "_", "")] = ""

                compound = cls(
                    smiles=compound_data.get(
                        "smiles", ""
                    ),  # blank default in case of solvents/reagents
                    idx=number_of_compounds,
                    reaction_smiles=reaction_table_dict.get("reaction_smiles", ""),
                    workbook=workbook,
                    demo="no",
                    reaction_component=component_type.capitalize(),
                    reaction_component_idx=len(component_list),
                    reload=True,
                )

                compound.compound_data.update(compound_data)
                if component_type == "solvent":
                    compound.add_solvent_sustainability_flags()

                component_lists[component_type].append(compound)

        return component_lists, units

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


@reaction_table_bp.route("/autoupdate_reaction_table", methods=["GET", "POST"])
# @workbook_member_required
def autoupdate_reaction_table():
    """
    I guess the idea is to take any smiles from the front end and generate the reaction table

    Needs to get reaction from db too, but ignore for now
    """
    # get user workbook
    demo = request.json.get("demo")
    tutorial = request.json.get("tutorial")
    workbook = None
    reaction = None
    polymer_indices = []
    if demo != "demo" and tutorial != "yes":
        workgroup = request.json.get("workgroup")
        workbook_name = request.json.get("workbook")
        workbook = services.workbook.get_workbook_from_group_book_name_combination(
            workgroup, workbook_name
        )
        reaction_id = request.json.get("reaction_id")
        reaction = services.reaction.get_from_reaction_id_and_workbook_id(
            reaction_id, workbook.id
        )

        polymer_indices = request.json.get("polymer_indices")

    default_units = {
        "limiting_reactant_table_number": 1,
        "main_product": -1,
        "mass_units": "mg",
        "polymerisation_type": "NAN",
        "amount_units": "mmol",
        "volume_units": "mL",
        "solvent_volume_units": "mL",
        "product_mass_units": "mg",
        "product_amount_units": "mmol",
    }

    reaction_smiles = request.json.get("reaction_smiles")
    if not reaction_smiles:
        return jsonify({"error": "Missing data!"})

    reactants_smiles_list, product_smiles_list = get_reactants_and_products_list(
        reaction_smiles
    )
    reactants = [
        SketcherCompound(
            smiles=x,
            idx=idx + 1,
            polymer_indices=polymer_indices,
            workbook=workbook,
            demo=demo,
            reaction_component="Reactant",
            reaction_component_idx=idx + 1,
        )
        for idx, x in enumerate(reactants_smiles_list)
    ]
    number_of_reactants = len(reactants)

    # add for reagent support
    # reagents = [
    #     SketcherCompound(
    #         smiles=x,
    #         idx=idx + 1 + number_of_reactants,
    #         polymer_indices=polymer_indices,
    #         workbook=workbook,
    #         demo=demo,
    #         reaction_component="Reagent"
    #     )
    #     for idx, x in enumerate(reactants_smiles_list)
    # ]
    # number_of_reagents = len(reactants)

    products = [
        SketcherCompound(
            smiles=y,
            idx=idx + 1 + number_of_reactants,
            polymer_indices=polymer_indices,
            workbook=workbook,
            demo=demo,
            reaction_component="Product",
            reaction_component_idx=idx + 1,
        )
        for idx, y in enumerate(product_smiles_list)
    ]
    number_of_products = len(products)

    # check errors first add reagents here for reagent support
    for compound_group in (reactants, products):
        # now handles co polymer, dummy atom and invalid molecule errors
        error = check_compound_errors(compound_group)
        if error:
            return error

        novel_compound_row = check_novel_compounds(compound_group)
        if novel_compound_row:
            # this should be changed i think
            return novel_compound_row

    identifiers = []

    # Solvents - keep solvents that are not novel compounds or are novel compounds within the current workbook
    if demo == "demo":
        sol_rows = services.solvent.get_default_list()
    else:
        sol_rows = services.solvent.get_workbook_list(workbook)

    r_class = None

    if not polymer_indices:
        polymer_indices = list()
        r_class = services.reaction_classification.classify_reaction(
            reactants_smiles_list, product_smiles_list
        )

    reaction_table_html = "reactions/_reaction_table.html"

    # Now it renders the reaction table template
    reaction_table = render_template(
        reaction_table_html,
        reactants=reactants,
        # reagents=reagents,
        number_of_reactants=number_of_reactants,
        number_of_products=number_of_products,
        number_of_reagents=0,
        identifiers=identifiers,
        reactant_table_numbers=[],
        products=products,
        units=default_units,
        # product_intended_dps=product_data["intended_dps"],
        reagent_table_numbers=[],
        reaction_table_data="",
        summary_table_data="",
        sol_rows=sol_rows,
        reaction=reaction,
        reaction_class=r_class,
        polymer_indices=polymer_indices,
    )
    return jsonify({"reactionTable": reaction_table})


@reaction_table_bp.route("/reload_reaction_table", methods=["GET", "POST"])
def reload_reaction_table():
    workbook = request.json.get("workbook")
    workgroup = request.json.get("workgroup")
    reaction_id = request.json.get("reaction_id")
    demo = request.json.get("demo")
    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup, workbook
    )
    reaction = services.reaction.get_from_reaction_id_and_workbook_id(
        reaction_id, workbook.id
    )
    compounds, units = SketcherCompound.from_reaction_table_dict(
        json.loads(reaction.reaction_table_data), workbook
    )

    if demo == "demo":
        sol_rows = services.solvent.get_default_list()
    else:
        sol_rows = services.solvent.get_workbook_list(workbook)

    # Now it renders the reaction table template
    reaction_table = render_template(
        "reactions/_reaction_table.html",
        reactants=compounds["reactant"],
        reagents=compounds["reagent"],
        solvents=compounds["solvent"],
        number_of_reactants=len(compounds["reactant"]),
        number_of_products=len(compounds["product"]),
        number_of_reagents=len(compounds["reagent"]),
        number_of_solvents=len(compounds["solvent"]),
        units=units,
        identifiers=[],
        reactant_table_numbers=[],
        products=compounds["product"],
        # product_intended_dps=product_data["intended_dps"],
        reagent_table_numbers=[],
        reaction_table_data="",
        summary_table_data="",
        sol_rows=sol_rows,
        reaction=reaction,
        reaction_class=reaction.reaction_class,
        reaction_classes=[],
        polymer_indices={},
    )
    return jsonify({"reactionTable": reaction_table})


def get_reactants_and_products_list(
    reaction_smiles: str,
) -> Tuple[List[str], List[str]]:
    """
    Process reactants and products strings to obtain lists of reactant and product SMILES.
    Converts words into SMILES symbols and identifies the format as either CXSMILES or SMILES.
    If ions are present, it processes the ionic SMILES accordingly.

    Args:
        reaction_smiles (str), the smiles of the reaction from the front end

    Returns:
        Tuple[List[str], List[str]]: A tuple containing lists of reactant, reagent and product SMILES.

    """
    # reactants_smiles = smiles_symbols(reactants)
    # products_smiles = smiles_symbols(products)

    # # form the reaction_smiles. Just from reactands and products to exclude reagents/other data
    # reaction_smiles = (reactants_smiles + ">>" + products_smiles).replace(",", ".")

    # we get the smiles straight from the sketcher to see if we have CXSmiles
    # sketcher_smiles = smiles_symbols(request.json.get("reaction_smiles"))

    # [OH-].[Na+]>>[Cl-].[Cl-].[Zn++] |f:0.1,2.3.4|
    if "|f:" in reaction_smiles:
        reaction_smiles += re.search(r" \|[^\|]*\|$", reaction_smiles).group()
        (
            reactants_smiles_list,
            products_smiles_list,
        ) = services.ions.reactants_and_products_from_ionic_cx_smiles(reaction_smiles)
    elif re.findall(r"(?<!\{)\+", reaction_smiles) or re.findall(
        r"(?<!\{)-", reaction_smiles
    ):
        # find "+" or "-" in reaction_smiles but not {-} or {+n} for polymers
        # CHANGE THIS FUNCTION TO HANDLE REAGENTS
        (
            reactants_smiles_list,
            products_smiles_list,
        ) = services.ions.reactants_and_products_from_ionic_smiles(reaction_smiles)
    # reactions with no ions - make rxn object directly from string
    else:
        # Split and pad to ensure two parts
        # change to 3 for reagent support and change split type >> to >
        parts = reaction_smiles.split(">>") + [""] * (
            2 - len(reaction_smiles.split(">>"))
        )

        # Process each part (add reagent_smiles_list_here)
        reactants_smiles_list, products_smiles_list = [
            x.split(".") if x else [] for x in parts
        ]

    return reactants_smiles_list, products_smiles_list


def get_compound_all_tables(smiles, workbook, polymer, demo):
    """
    Retrieves a compound from the database, checking the Compound, NovelCompound, and PolymerNovelCompound tables.

    Args:
        smiles (str or list): The SMILES string for the compound.
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

    compound_data["molecular_weights"] = []
    compound_data["names"] = []
    compound_data["hazards"] = []
    compound_data["densities"] = []
    compound_data["primary_keys"] = []

    if isinstance(compound, models.PolymerNovelCompound):
        molecular_weight = services.polymer_novel_compound.get_repeat_unit_weights(
            compound.id, compound.workbook
        )
    else:
        molecular_weight = (
            float(compound.molec_weight) if compound.molec_weight != "" else 0
        )

    compound_data["molecular_weights"].append(molecular_weight)

    compound_name = compound.name if compound.name != "" else "Not found"
    compound_data["names"].append(compound_name)

    compound_hazard = (
        compound.hphrase if compound.hphrase != "No hazard codes found" else "Unknown"
    )
    compound_data["hazards"].append(compound_hazard)

    compound_density = compound.density if compound.density != "" else "-"
    compound_data["densities"].append(compound_density)

    if novel_compound:
        compound_data["primary_keys"].append((compound.name, compound.workbook))
    else:
        compound_data["primary_keys"].append(compound.id)


@reaction_table_bp.route("/_save_reaction_note", methods=["POST"])
@login_required
@reaction_table_bp.doc(security="sessionAuth")
def save_reaction_note():
    """
    Saves an reaction_note to the reaction object

    Returns:
        flask.Response: A JSON response with the reaction_note object
    """
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
