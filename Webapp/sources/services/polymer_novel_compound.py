import re
from datetime import datetime
from typing import Optional, Tuple

from flask import abort, jsonify, request
from flask_login import current_user
from psmiles import PolymerSmiles
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
        .filter(func.lower(models.PolymerNovelCompound.name) == primary_key[0].lower())
        .join(models.WorkBook)
        .filter(models.WorkBook.id == primary_key[1])
        .first()
    )[0]


def from_name_and_workbook(
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


def from_smiles_and_workbook(
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

        check_in_workbook = services.polymer_novel_compound.from_name_and_workbook(
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
        check_in_workbook = services.polymer_novel_compound.from_smiles_and_workbook(
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


def extract_inside_brackets(
    s: str,
    start_index: int,
) -> tuple[str, int]:
    """
    Extracts contents of brackets to identify branching in polymer SMILES.
    Handles nested branching.

    Args:
        s - whole SMILES string
        start_index - index of first bracket

    Returns:
        The branch SMILES string and end index, or -1 if no branch found
    """
    balance = 0
    inner_content = ""
    for i in range(start_index, len(s)):
        if s[i] == "(":
            balance += 1
        elif s[i] == ")":
            balance -= 1

        inner_content += s[i]  # Capture the character

        if balance == 0:  # Found matching closing parenthesis
            return inner_content, i

    return inner_content, -1  # Return -1 if no matching parenthesis found


def check_bracket_balance(s):
    balance = 0
    for i in range(len(s)):
        if s[i] == "(":
            balance += 1
        elif s[i] == ")":
            balance -= 1

    return balance


def find_polymer_repeat_unit(
    compound: str,
) -> str:
    """
    Identify the repeat group within a polymer SMILES string.
    Must be done before smiles is cleaned.
    Uses {-} and {+n} to find repeat unit.
    """
    start_marker = compound.find("{-}")
    end_marker = compound.find("{+n}")

    # if either {-} or {+} is not found
    if start_marker == -1 or end_marker == -1:
        return ""

    # START is the letter before {-}
    if compound[start_marker - 1].isupper():
        result = "*" + compound[start_marker - 1] + compound[start_marker + 3 :]
    elif (
        compound[start_marker - 1] == "]"
    ):  # for polymers with [] groups e.g C[SiH2]{-}CC{+n}C
        i = 2
        while compound[start_marker - i] != "[" and start_marker - i > 0:
            i += 1  # search backwards
        result = "*" + compound[start_marker - i :].replace("{-}", "")
    else:  # for polymers with rings e.g *C1{-}CC{+n}(CCC1)*
        i = 2
        while not compound[start_marker - i].isupper() and start_marker - i > 0:
            i += 1  # search backwards
        result = "*" + compound[start_marker - i :].replace("{-}", "")

    # END is the character before {+n}
    if ")" not in result[end_marker:]:  # no branching at end of SRU
        return result[: end_marker - 3] + "*"

    if result[end_marker + 1] == "(":  # last atom of SRU has a branch
        result = result[: end_marker - 3]
        branch, close_paren = extract_inside_brackets(compound, end_marker + 4)
        if close_paren != -1:
            result += branch

        # check for more branching on last atom of SRU
        while close_paren + 1 < len(compound) and compound[close_paren + 1] == "(":
            branch, close_paren = extract_inside_brackets(compound, close_paren + 1)
            if close_paren != -1:
                result += branch
    else:
        result = result[: result.find("{+n}")]

    # if end is inside brackets, ----------------------
    if check_bracket_balance(compound[:end_marker]) == 0:
        return result + "*"

    section = []
    parts = compound[end_marker + 4 :].split(")")
    for i, part in enumerate(parts):
        if "(" in part and i < len(parts) - 1:
            parts[i] = ")" + part + ")" + parts[i + 1]
            parts.pop(i + 1)
            section.append(parts[i])
        else:
            section.append(")" + part)
    idx = check_bracket_balance(compound[:end_marker])

    result = result + "*" + "".join(section[-idx:])

    return result


def clean_polymer_smiles(
    compound: str,
) -> str:
    """
    Removes polymer symbols from the string generated by sketchers. last one removes dangling bonds
    """
    return compound.replace("{+n}", "").replace("{-}", "").replace("-", "")


def canonicalise(smiles: str) -> str:
    """
    Canonicalise a polymer SMILES string. Not for use when endgroups are known (must have two '*'s)
    """
    try:
        ps = PolymerSmiles(smiles)
        smiles = ps.canonicalize
        return str(smiles).replace("[*]", "*")  # change dangling bonds label
    except Exception as e:
        return f"Error: {str(e)}"


def find_canonical_repeat(smiles: str) -> str:
    """
    Find repeat unit, and canonicalise a polymer smiles string.
    REMOVES END GROUPS
    Args:
        smiles - the unedited smiles from ketcher output. e.g. CC{-}C{+n}C
    Returns:
        canon_smiles - the canonicalised repeat unit. e.g *C*
    """
    if "[*]" in smiles:  # block dummy atoms like R groups
        return "dummy"
    smiles = find_polymer_repeat_unit(smiles)

    smiles = smiles.replace("/", "").replace("\\", "")  # ignore stereochemistry

    smiles = canonicalise(smiles)

    return smiles
