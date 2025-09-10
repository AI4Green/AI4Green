from typing import Dict, List, Optional, Union
from urllib.parse import quote
from urllib.request import urlopen

from flask import jsonify, render_template
from rdkit import Chem
from sources import models, services
from sources.extensions import db
from sqlalchemy import func

"""
For all things associated with the Compound table in the database
"""


def count() -> int:
    """
    Gets the number of compounds in the database

    Returns:
        The number of compounds in the database
    """
    return db.session.query(models.Compound).count()


def get_compound_data_error_reports() -> List[models.CompoundDataErrorReport]:
    """
    Gets a list of compound database error reports in the database

    Returns:
         List of all compound database error reports
    """
    # Add time cut off if list grows too large
    return db.session.query(models.CompoundDataErrorReport).all()


def get_smiles(primary_key: int) -> str:
    """
    Gets the smiles from the compound primary key
    """
    return (
        db.session.query(models.Compound.smiles)
        .filter(models.Compound.id == primary_key)
        .first()
    )[0]


def get(primary_key: int) -> models.Compound:
    """Returns the Compound record that matches the primary key ID"""
    return (
        db.session.query(models.Compound)
        .filter(models.Compound.id == primary_key)
        .first()
    )


def from_smiles(smiles: str) -> Optional[models.Compound]:
    """
    Retrieve a compound from the database by converting SMILES to InChI.

    Args:
        smiles (str): The SMILES representation of the compound.

    Returns:
        models.Compound: A Compound object corresponding to the provided SMILES.
    """
    inchi = services.all_compounds.smiles_to_inchi(smiles)
    if not inchi:
        return None
    return from_inchi(inchi)


def from_inchi(inchi: str) -> models.Compound:
    """
    Retrieve a compound from the database based on its InChI.

    Args:
        inchi (str): The InChI (International Chemical Identifier) of the compound.

    Returns:
        models.Compound: The compound object retrieved from the database.
                        Returns None if no matching compound is found.
    """
    return (
        db.session.query(models.Compound).filter(models.Compound.inchi == inchi).first()
    )


def from_cas(cas: str) -> models.Compound:
    """
    Retrieve a compound from the database based on its CAS number

    Args:
        CAS - the chemical abstract service registry number of a compound

    Returns:
        models.Compound: The compound object retrieved from the database.
                        Returns None if no matching compound is found.
    """
    return db.session.query(models.Compound).filter(models.Compound.cas == cas).first()


def from_name(name: str) -> models.Compound:
    """
    Retrieve a compound from the database based on its name

    Args:
        name - the name of a compound

    Returns:
        models.Compound: The compound object retrieved from the database.
                        Returns None if no matching compound is found.
    """
    return (
        db.session.query(models.Compound)
        .filter(func.lower(models.Compound.name) == name.lower())
        .first()
    )


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
            self.is_novel_compound = True

        else:
            self.is_novel_compound = novel_compound
            self.inchi = compound.inchi
            get_compound_data(self.compound_data, compound, novel_compound)

    def check_for_polymer(self, polymer_indices, reaction_smiles):
        if polymer_indices is not None:
            self.check_polymer_indices_for_polymer(polymer_indices)
        # else:
        #     self.check_reaction_smiles_for_polymer(reaction_smiles)

    def check_polymer_indices_for_polymer(self, polymer_indices):
        """
        Set is_polymer to True if idx in polymer indices and process smiles
        """
        if self.idx in polymer_indices:
            self.is_polymer = True
            self.smiles = services.polymer_novel_compound.find_canonical_repeats(
                self.smiles
            )

    def check_reaction_smiles_for_polymer(self, reaction_smiles):
        (
            reactant_smiles,
            product_smiles,
        ) = services.reaction_table.get_reactants_and_products_list(reaction_smiles)
        smiles_list = (
            reactant_smiles if self.reaction_component == "Reactant" else product_smiles
        )
        # only check reactant and products for polymers
        if self.reaction_component in ["Reactant", "Product"]:
            print(smiles_list, self.reaction_component_idx)
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
            "reactions/_novel_compound.html",
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
                    try:
                        compound_data[key.replace(component_type + "_", "")] = value[
                            idx
                        ]
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

        return component_lists

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
