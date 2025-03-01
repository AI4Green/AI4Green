import re
from typing import List, Optional, Tuple

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
        db.session.query(models.NovelCompound.smiles)
        .filter(func.lower(models.NovelCompound.name) == primary_key[0].lower())
        .join(models.WorkBook)
        .filter(models.WorkBook.id == primary_key[1])
        .first()
    )[0]


def from_name_and_workbook(
    name: str, workbook: models.WorkBook
) -> models.NovelCompound:
    """
    Retrieves a novel compound by name and workbook.

    Args:
        name: Compound name.
        workbook: Workbook model.

    Returns:
        NovelCompound model.
    """
    return (
        db.session.query(models.NovelCompound)
        .filter(func.lower(models.NovelCompound.name) == name.lower())
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook.id)
        .first()
    )


def from_inchi_and_workbook(
    inchi: str, workbook: models.WorkBook
) -> models.NovelCompound:
    """
    Retrieves a novel compound by name and workbook.

    Args:
        inchi: Compound inchi.
        workbook: Workbook model.

    Returns:
        NovelCompound model.
    """
    return (
        db.session.query(models.NovelCompound)
        .filter(models.NovelCompound.inchi == inchi)
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook.id)
        .first()
    )


def from_cas_and_workbook(cas: str, workbook: models.WorkBook) -> models.NovelCompound:
    """
    Retrieves a novel compound by CAS and workbook.

    Args:
        cas - the chemical abstract service registry number of a novel compound
        workbook - Workbook model.

    Returns:
        NovelCompound model.
    """
    return (
        db.session.query(models.NovelCompound)
        .filter(models.NovelCompound.cas == cas)
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook.id)
        .first()
    )


def add(
    name: str,
    cas: Optional[str],
    mol_formula: str,
    mol_weight: float,
    density: float,
    concentration: float,
    hazards: str,
    smiles: str,
    inchi: Optional[str],
    inchi_key: Optional[str],
    workbook_id: int,
) -> models.NovelCompound:
    """
    Creates a novel compound in the database.

    Args:
        name: Compound name.
        cas: CAS number.
        mol_formula: Molecule formula.
        mol_weight: Molecular weight.
        density: Density.
        concentration: Concentration.
        hazards: Hazard codes.
        smiles: SMILES notation.
        inchi: InChI notation.
        inchi_key: InChIKey.
        workbook_id: Workbook ID.

    Returns:
        NovelCompound model.
    """
    nc = models.NovelCompound(
        name=name,
        cas=cas,
        molec_formula=mol_formula,
        molec_weight=mol_weight,
        density=density,
        concentration=concentration,
        hphrase=hazards,
        smiles=smiles,
        inchi=inchi,
        inchikey=inchi_key,
        workbook=workbook_id,
    )
    db.session.add(nc)
    db.session.commit()
    return nc


def reform_novel_compound_primary_key(primary_key: str) -> Tuple:
    """
    Converts a novel primary key to a tuple from the string returned from the frontend HTML
    Used for polymers too

    Args:
        primary_key - the primary key as a string. e.g., ('pixie dust', 1)

    Returns:
        A tuple of (compound_name, workbook_id)
    """
    if len(primary_key) > 350:
        abort(
            413
        )  # content too large. Exceeds max workbook name length + max novel compound name length
    # extract the name from between either the single quotes or double quotes.
    compound_name = extract_compound_name(primary_key)
    # extract the workbook id from after the end of the single or double quote then comma
    workbook_id = int(re.search(r"""["'], (\d+)""", primary_key).group(1))
    return compound_name, workbook_id


def extract_compound_name(primary_key: str) -> Optional[str]:
    """
    Extracts and unescapes the compound name from the primary_key string.
    Uses two regex patterns to handle both single and double quotes.
    example formats - : "('compound name', 1)" or '("compound's name", 1)'
    Args:
        primary_key (str): The input string containing the compound name.
    Returns:
        str or None: The extracted compound name, or None if no match is found.
    """
    # One regex for single quotes and one for double quotes
    pattern_single = r"\('((?:\\'|[^'])*)', \d"
    pattern_double = r'\("((?:\\"|[^"])*)", \d'
    # Try matching with single quotes pattern
    name_match = re.search(pattern_single, primary_key)
    if name_match:
        compound_name = name_match.group(1)
        return compound_name

    # Try matching with double quotes pattern if single quotes pattern fails
    name_match = re.search(pattern_double, primary_key)
    if name_match:
        compound_name = name_match.group(1)
        return compound_name

    # Return None if no match is found
    return None


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
        self.smiles = sanitise_user_input(request.form["smiles"])
        self.cas = sanitise_user_input(request.form["cas"])
        self.mol_weight = self._get_mol_weight()
        self.feedback = ""
        self.validation = ""
        self.validation_functions = [
            self.validate_name_is_allowed,
            self.validate_name_is_unique,
            self.validate_cas_is_allowed,
            self.validate_cas_is_unique,
            self.validate_hazards,
            self.validate_numerical_properties,
            self.validate_smiles,
            self.validate_structure_is_unique,
        ]

        # if SMILES is provided mol_formula, inchi, and inchi_key are calculated
        self.mol_formula = ""
        self.inchi = None
        self.inchi_key = None
        self.calculate_molecule_identifiers()

    def _get_mol_weight(self) -> Optional[float]:
        """Gets molecular weight from the request form or calculate from smiles if not provided"""
        if request.form["molWeight"] and request.form["molWeight"] != "":
            return float(request.form["molWeight"])
        elif request.form["smiles"] and Chem.MolFromSmiles(self.smiles):
            return services.all_compounds.mol_weight_from_smiles(self.smiles)
        else:
            return None

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
        if not self.name:
            self.feedback = "Compound requires a name"
            self.validation = "failed"

        # maximum name length
        elif len(self.name) > 200:
            self.feedback = "Name must be under 200 characters long"
            self.validation = "failed"

        # only certain symbols allowed
        if not re.fullmatch(r"[A-Za-z0-9α-ωΑ-Ω+\-.,/(){}\[\]·'\s]+", self.name):
            invalid_symbols = self.find_invalid_symbols(self.name)
            self.feedback = (
                f"Name contains invalid symbols {', '.join(invalid_symbols)}"
            )
            self.validation = "failed"

    @staticmethod
    def find_invalid_symbols(chemical_name: str):
        """
        Identify the invalid symbols in the provided chemical name.

        Args:
             chemical_name: The chemical name we are finding invalid symbols in.
        Returns:
            A set of invalid symbols found in the chemical name.
        """
        valid_symbols = set(
            "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789αβγδεζηθικλμνξοπρστυφχψωΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΠΡΣΤΥΦΧΨΩ+-.,/(){}[]·' "
        )

        invalid_symbols = set()
        for char in chemical_name:
            if char not in valid_symbols:
                invalid_symbols.add(char)
        return invalid_symbols

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

        check_in_workbook = services.novel_compound.from_name_and_workbook(
            self.name, self.workbook
        )
        if check_in_workbook:
            self.feedback = (
                "There is already a compound with this name in this workbook"
            )
            self.validation = "failed"

    def validate_cas_is_allowed(self):
        """
        Validates the CAS number from the request. Checks it fits the expected CAS pattern with a regular expression
        """
        cas = sanitise_user_input(request.form["cas"])
        if cas:
            cas_regex = r"^[0-9]{1,7}-\d{2}-\d$"
            if not re.findall(cas_regex, cas):
                self.feedback = "Invalid CAS number"
                self.validation = "failed"

    def validate_cas_is_unique(self):
        """
        Checks if the compound CAS number is unique within the compound database and the given workbook.
        """
        if not self.cas:
            return
        check_in_compound_database = services.compound.from_cas(self.cas)
        if check_in_compound_database:
            self.feedback = (
                "A compound with this CAS number is already in the compound database"
            )
            self.validation = "failed"

        check_in_workbook = services.novel_compound.from_cas_and_workbook(
            self.cas, self.workbook
        )
        if check_in_workbook:
            self.feedback = (
                "A compound with this CAS number is already in this workbook"
            )
            self.validation = "failed"

    def validate_structure_is_unique(self):
        """
        Optional validation if the user supplies a SMILES string or draws the compound.
        Confirms the structure is not already in the database in either the user's workbook or the compound database.
        """
        if not self.inchi:
            return
        check_in_compound_db = services.compound.from_inchi(self.inchi)
        if check_in_compound_db:
            self.feedback = f"There is already a compound with this structure in the Compound database with the name: {check_in_compound_db.name}"
            self.validation = "failed"

        check_in_workbook = services.novel_compound.from_inchi_and_workbook(
            self.inchi, self.workbook
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
            entry is None or services.utils.check_positive_number(entry)
            for entry in numerical_values
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
        self.inchi = Chem.MolToInchi(mol)
        self.inchi_key = Chem.MolToInchiKey(mol)

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


def all_from_workbook(workbook: models.WorkBook) -> List[models.NovelCompound]:
    """
    Retrieves all novel compounds in a workbook.

    Args:
        workbook: Workbook model.

    Returns:
        A list of novel compounds.
    """
    return (
        db.session.query(models.NovelCompound)
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook.id)
        .all()
    )
