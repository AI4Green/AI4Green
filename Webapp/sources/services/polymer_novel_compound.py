from datetime import datetime
from typing import Optional, Tuple

from flask import abort, request
from flask_login import current_user
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from sources import models, services
from sources.auxiliary import sanitise_user_input
from sources.extensions import db
from sqlalchemy import func


def get_smiles(primary_key: Tuple[str, int], person: models.Person = None) -> str:
    """
    Gets the novel compound's SMILES string from the primary key if the entry has the SMILES attribute
    Args:
        primary_key - the primary key is a tuple in the format [compound_name, workgroup.id]
        person - the person we are checking for access rights to the novel compound

    Returns:
        The SMILES string corresponding to the primary key or None
    """
    primary_key = (primary_key[0], int(primary_key[1]))
    workbook = services.workbook.get(primary_key[1])
    if (person and person not in workbook.users) or (
        not person and current_user.Person not in workbook.users
    ):
        abort(401)

    return (
        db.session.query(models.PolymerNovelCompound.smiles)
        .filter(models.PolymerNovelCompound.name == primary_key[0].lower())
        .join(models.WorkBook)
        .filter(models.WorkBook.id == primary_key[1])
        .first()
    )[0]


def get_polymer_novel_compound_from_name_and_workbook(
    name: str, workbook: models.WorkBook
) -> models.PolymerNovelCompound:
    """
    Retrieves a novel compound by name and workbook.

    Args:
        name: Compound name.
        workbook: Workbook model.

    Returns:
        PolymerNovelCompound model.
    """
    return (
        db.session.query(models.PolymerNovelCompound)
        .filter(func.lower(models.PolymerNovelCompound.name) == name.lower())
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook.id)
        .first()
    )


def get_polymer_novel_compound_from_smiles_and_workbook(
    smiles: str, workbook: models.WorkBook
) -> models.PolymerNovelCompound:
    """
    Retrieves a novel compound by name and workbook.

    Args:
        smiles: Compound smiles.
        workbook: Workbook model.

    Returns:
        PolymerNovelCompound model.
    """
    return (
        db.session.query(models.PolymerNovelCompound)
        .filter(models.PolymerNovelCompound.smiles == smiles)
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook.id)
        .first()
    )


def add(
    name: str,
    mol_formula: str,
    mol_weight: float,
    density: float,
    concentration: float,
    hazards: str,
    smiles: str,
    workbook_id: int,
) -> models.PolymerNovelCompound:
    """
    Creates a novel compound in the database.

    Args:
        name: Compound name.
        mol_formula: Molecule formula.
        mol_weight: Molecular weight.
        density: Density.
        concentration: Concentration.
        hazards: Hazard codes.
        smiles: SMILES notation.
        workbook_id: Workbook ID.
        current_time: Current timestamp.

    Returns:
        PolymerNovelCompound model.
    """
    nc = models.PolymerNovelCompound(
        name=name,
        molec_formula=mol_formula,
        molec_weight=mol_weight,
        density=density,
        concentration=concentration,
        hphrase=hazards,
        smiles=smiles,
        workbook=workbook_id,
    )
    db.session.add(nc)
    db.session.commit()
    return nc


class NewNovelCompound:
    def __init__(self, workbook):
        self.workbook = workbook
        self.name = sanitise_user_input(request.form["name"])
        self.density = (
            float(request.form["density"]) if request.form["density"] != "" else None
        )
        self.concentration = (
            float(request.form["concentration"])
            if request.form["concentration"] != ""
            else None
        )
        self.mol_weight = (
            float(request.form["molWeight"])
            if request.form["molWeight"] != ""
            else None
        )
        self.hazard_codes = sanitise_user_input(request.form["hPhrase"])
        self.smiles = sanitise_user_input(request.form["smiles"])  # canonicalised

        self.feedback = ""
        self.validation = ""
        self.validation_functions = [
            self.validate_name_is_allowed,
            self.validate_name_is_unique,
            self.validate_hazards,
            self.validate_numerical_properties,
            self.validate_smiles,
            self.validate_structure_is_unique,
        ]

        # if SMILES is provided mol_formula is calculated
        self.mol_formula = ""
        self.calculate_molecule_identifiers()

    def validate(self):
        for validation_function in self.validation_functions:
            validation_function()
            if self.validation == "failed":
                print("New compound validation failed due to: ", self.feedback)
                return

        self.feedback = "Compound added to the database"
        return "validation_successful"

    def validate_name_is_allowed(self):
        """
        Validates and sanitizes the compound name from the request.
        """
        if len(self.name) > 200:
            self.feedback = "Name must be under 200 characters long"
            self.validation = "failed"
        if not self.name:
            self.feedback = "Compound requires a name"
            self.validation = "failed"

    def validate_name_is_unique(self):
        """
        Checks if the compound name is unique within the compound database and the given workbook.
        """
        check_in_compound_database = services.compound.from_name(self.name)
        if check_in_compound_database:
            self.feedback = (
                "There is already a compound with this name in the compound database"
            )
            self.validation = "failed"

        check_in_workbook = services.polymer_novel_compound.get_polymer_novel_compound_from_name_and_workbook(
            self.name, self.workbook
        )
        if check_in_workbook:
            self.feedback = (
                "There is already a compound with this name in this workbook"
            )
            self.validation = "failed"

    def validate_structure_is_unique(self):
        """
        Optional validation if the user supplies a SMILES string or draws the compound.
        Confirms the structure is not already in the database in the user's workbook.
        Uses SMILES not inchi.
        """
        if not self.smiles:
            return
        check_in_workbook = services.polymer_novel_compound.get_polymer_novel_compound_from_smiles_and_workbook(
            self.smiles, self.workbook
        )
        if check_in_workbook:
            self.feedback = f"There is already a compound with this structure saved in this workbook with the name: {check_in_workbook.name}"
            self.validation = "failed"

    def validate_numerical_properties(self):
        """
        Density, concentration, and molecular weight are optional properties that if provided must be a positive number.
        This function validates this data or fails the validation and provides feedback.
        """
        numerical_values = [self.density, self.concentration, self.mol_weight]
        if all(
            entry is None or check_positive_number(entry) for entry in numerical_values
        ):
            return
        self.feedback = "Molecular weight, density, and concentration must be empty or a positive number"
        self.validation = "failed"

    def validate_smiles(self):
        """
        If the SMILES string is provided confirm that it can generate a structure in rdkit.
        Validation ignored for compounds drawn in the sketcher to prevent blocking users if there is issue with rdkit.
        """
        if self.smiles and request.form["source"] != "sketcher":
            mol = Chem.MolFromSmiles(self.smiles)
            if mol is None:
                self.feedback = (
                    f"Failed to generate molecule from SMILES string: {self.smiles}"
                )
                self.validation = "failed"

    def calculate_molecule_identifiers(self):
        """
        Calculate additional molecule identifiers if SMILES is present.
        """
        mol = Chem.MolFromSmiles(self.smiles)
        if mol is None or mol.GetNumAtoms() == 0:
            return
        self.mol_formula = rdMolDescriptors.CalcMolFormula(mol)

    def validate_hazards(self):
        """
        Check hazard codes are in the correct format.
        """
        if not self.hazard_codes:
            self.hazard_codes = "Unknown"
            return
        hazard_codes_ls = self.hazard_codes.split("-")
        for hazard_code in hazard_codes_ls:
            hazard_match = services.hazard_code.get(hazard_code)
            if hazard_match is None:
                self.feedback = f'Hazard code "{hazard_code}" is invalid. Must be valid hazard code and formatted correctly. e.g., H200-H301.'
                self.validation = "failed"


def check_positive_number(s: float) -> bool:
    """Checks the entry is a positive number."""
    try:
        return s >= 0
    except (ValueError, TypeError):
        return False
