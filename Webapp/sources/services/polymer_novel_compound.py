import re
from typing import Tuple

from flask import abort, json, request
from flask_login import current_user
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from sources import models, services
from sources.auxiliary import sanitise_user_input
from sources.extensions import db
from sqlalchemy import exists, func


def get_smiles(primary_key: Tuple[str, int], person: models.Person = None) -> list[str]:
    """
    Gets the repeat units SMILES strings from the PolymerNovelCompound primary key by searching in the PolymerRepeatUnit table
    Args:
        primary_key - the primary key is a tuple in the format [compound_name, workgroup.id]
        person - the person we are checking for access rights to the novel compound

    Returns:
        A list of the SMILES strings corresponding to the primary key or None
    """
    primary_key = (primary_key[0], int(primary_key[1]))
    workbook = services.workbook.get(primary_key[1])
    if (person and person not in workbook.users) or (
        not person and current_user.Person not in workbook.users
    ):
        abort(401)

    query = (
        db.session.query(models.PolymerRepeatUnit.smiles)
        .filter(func.lower(models.PolymerRepeatUnit.name) == primary_key[0].lower())
        .join(models.WorkBook)
        .filter(models.WorkBook.id == primary_key[1])
        .all()
    )
    return [q[0] for q in query]


def get_repeat_unit_weights(polymer_id: int, workbook: int) -> list:
    """
    Retrieves a list of molec weights of the repeat units, from the polymer id.

    Args:
        polymer_id: PolymerNovelCompound.id
        workbook: PolymerNovelCompound.workbook

    Returns:
        PolymerNovelCompound model.
    """
    query = (
        db.session.query(models.PolymerRepeatUnit.molec_weight)
        .filter(models.PolymerRepeatUnit.polymer_id == polymer_id)
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook)
        .all()
    )
    return [q[0] for q in query]


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
    smiles: list, workbook: models.WorkBook
) -> models.PolymerNovelCompound:
    """
    Retrieves a novel compound by repeat units and workbook.

    Args:
        smiles: Repeat unit smiles as list.
        workbook: Workbook model.

    Returns:
        PolymerNovelCompound model.
    """
    query = (
        db.session.query(models.PolymerNovelCompound)
        .filter(
            # Ensure no repeat units not in smiles list
            ~exists().where(
                (models.PolymerRepeatUnit.polymer_id == models.PolymerNovelCompound.id)
                & (~models.PolymerRepeatUnit.smiles.in_(smiles))
            ),
            # Ensure all repeat units are from the smiles list and count matches
            db.session.query(func.count(models.PolymerRepeatUnit.id))
            .filter(
                models.PolymerRepeatUnit.polymer_id == models.PolymerNovelCompound.id
            )
            .scalar_subquery()
            == len(smiles),
        )
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook.id)
        .first()
    )

    return query


def add(
    name: str,
    mol_formula: list,
    mol_weight: list,
    density: float,
    concentration: float,
    hazards: str,
    smiles: list,
    workbook_id: int,
) -> models.PolymerNovelCompound:
    """
    Creates a novel compound in the database, and adds repeat units to database

    Args:
        name: Compound name.
        mol_formula: Molecule formula.
        mol_weight: Molecular weight.
        density: Density.
        concentration: Concentration.
        hazards: Hazard codes.
        smiles: list of strings in SMILES notation.
        workbook_id: Workbook ID.

    Returns:
        PolymerNovelCompound model.
    """
    nc = models.PolymerNovelCompound(
        name=name,
        density=density,
        concentration=concentration,
        hphrase=hazards,
        workbook=workbook_id,
    )
    db.session.add(nc)
    db.session.commit()

    polymer = from_name_and_workbook(name, services.workbook.get(workbook_id))

    for idx in range(len(smiles)):
        nr = models.PolymerRepeatUnit(
            name=name,
            polymer_id=polymer.id,
            polymer=polymer,
            smiles=smiles[idx],
            molec_weight=mol_weight[idx],
            molec_formula=mol_formula[idx],
            workbook=workbook_id,
        )
        db.session.add(nr)
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
            [float(mw) for mw in json.loads(request.form["molWeight"])]
            if request.form["molWeight"] != ""
            else None
        )
        self.hazard_codes = sanitise_user_input(request.form["hPhrase"])
        self.smiles = [
            sanitise_user_input(smiles) for smiles in json.loads(request.form["smiles"])
        ]  # canonicalised

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
        numerical_values = [self.density, self.concentration]
        numerical_values.extend(self.mol_weight)
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
        mol_formula = []
        for smiles in self.smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None or mol.GetNumAtoms() == 0:
                mol_formula.append("")  # TODO: is this gonna cause issues?
                continue
            mol_formula.append(rdMolDescriptors.CalcMolFormula(mol))
        self.mol_formula = mol_formula

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


def check_bracket_balance(smiles):
    """
    Returns the number of opening brackets compared to closing brackets in a SMILES string.
    """
    balance = 0
    for i in range(len(smiles)):
        if smiles[i] == "(":
            balance += 1
        elif smiles[i] == ")":
            balance -= 1

    return balance


def find_polymer_repeat_units(
    compound: str,
) -> list[str]:
    """
    Identify the repeat group within a polymer SMILES string.
    Must be done before smiles is cleaned.
    Uses {-} and {+n} to find repeat unit.
    """
    start_marker = [marker.start() for marker in re.finditer("{-}", compound)]
    end_marker = [marker.start() for marker in re.finditer(re.escape("{+n}"), compound)]

    # if either {-} or {+} is not found
    if start_marker == -1 or end_marker == -1:
        return [""]

    results = []
    for i in range(len(start_marker)):
        # START is the letter before {-}
        if compound[start_marker[i] - 1].isupper():
            result = (
                "*" + compound[start_marker[i] - 1] + compound[start_marker[i] + 3 :]
            )
        elif (
            compound[start_marker[i] - 1] == "]"
        ):  # for polymers with [] groups e.g C[SiH2]{-}CC{+n}C
            x = 2
            while compound[start_marker[i] - x] != "[" and start_marker[i] - x > 0:
                x += 1  # search backwards
            result = "*" + compound[start_marker[i] - x :].replace("{-}", "")
        else:  # for polymers with rings e.g *C1{-}CC{+n}(CCC1)*
            x = 2
            while (
                not compound[start_marker[i] - x].isupper() and start_marker[i] - x > 0
            ):
                x += 1  # search backwards
            result = "*" + compound[start_marker[i] - x :].replace("{-}", "")

        diff = len(compound) - len(result)

        # END is the character before {+n}
        if i + 1 < len(end_marker):
            if (
                ")" not in result[end_marker[i] - diff : end_marker[i + 1] - diff]
            ):  # no branching at end of SRU
                results.append(result[: end_marker[i] - diff] + "*")
                continue
        else:  # last repeat unit
            if ")" not in result[end_marker[i] - diff :]:  # no branching at end of SRU
                results.append(result[: end_marker[i] - diff] + "*")
                continue

        if result[end_marker[i] - diff + 4] == "(":  # last atom of SRU has a branch
            result = result[: end_marker[i] - diff]
            branch, close_paren = extract_inside_brackets(compound, end_marker[i] + 4)
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
        if check_bracket_balance(compound[start_marker[i] : end_marker[i]]) == 0:
            results.append(result + "*")
            continue

        section = []
        parts = compound[end_marker[i] + 4 :].split(")")
        for p, part in enumerate(parts):
            if "(" in part and i < len(parts) - 1:
                while check_bracket_balance(parts[p]) > 0:
                    parts[p] = ")" + parts[p] + ")" + parts[p + 1]
                    parts.pop(p + 1)
                section.append(parts[p])
            else:
                section.append(")" + part)

        idx = check_bracket_balance(compound[: end_marker[i]])

        result = result + "*" + "".join(section[-idx:])

        results.append(result)

    return results


def clean_polymer_smiles(
    compound: str,
) -> str:
    """
    Removes polymer symbols from the string generated by sketchers. last one removes dangling bonds
    """
    return compound.replace("{+n}", "").replace("{-}", "").replace("-", "")


def identify_symmetry(
    mol: Chem.rdchem.Mol,
    anchors: list,
) -> str:
    """
    identify symmetry in a molecule with two atoms and reduce

    :param mol: mol with two atoms in backbone e.g *CC*
    :param anchors: list of atom indices for the neighbours of dummy atoms
    :return: canonicalised and reduced smiles

    """
    # check for symmetry
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))

    if ranks[anchors[0]] != ranks[anchors[1]]:
        # no symmetry
        return Chem.MolToSmiles(mol)

    # break bond in symmetry plane and replace with dummy atoms. *CC* > *C*.*C*
    fragmented_mol = Chem.FragmentOnBonds(
        mol,
        [mol.GetBondBetweenAtoms(anchors[0], anchors[1]).GetIdx()],
        dummyLabels=[(0, 0)],
    )
    frags = Chem.GetMolFrags(fragmented_mol, asMols=True, sanitizeFrags=True)

    # return the first fragment
    return Chem.MolToSmiles(frags[0])


def all_ring_bonds(mol):
    """
    check if all bonds in mol are within a ring

    :param mol: mol
    :return: True if all bonds are within a ring
    """

    for b in mol.GetBonds():
        ring_bond = b.HasProp("ring")

        # Check if the bond is attached to a dummy atom
        attached_to_dummy = (
            b.GetBeginAtom().GetAtomicNum() == 0 or b.GetEndAtom().GetAtomicNum() == 0
        )

        if not (ring_bond or attached_to_dummy):
            # standard linear backbone
            return False

    # everything is in a ring
    return True


def get_repeating_substructure(mol):
    """
    identifies repeating substructures in a cyclic mol, based on rotational symmetry and returns substructure

    :param mol: cyclic mol object
    :return: reduced cyclic mol object, boolean indicating whether repeats were found
    """
    Chem.SanitizeMol(mol)

    # Get topological symmetry ranks
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    ssr = Chem.GetSymmSSSR(mol)

    for ring in ssr:
        ring_indices = sorted(list(ring))
        n = len(ring_indices)
        ring_ranks = [ranks[i] for i in ring_indices]

        # Find smallest repeating period, i, to divide length n perfectly
        period = n
        for i in range(1, n // 2 + 1):
            if n % i == 0:
                # Check if the ring ranks are a repetition of the first i elements
                if ring_ranks[:i] * (n // i) == ring_ranks:
                    period = i
                    break

        # if repeats found
        if period < n:
            unit_ring_indices = ring_indices[:period]

            # Capture ring segment + substituents
            final_atom_set = set(unit_ring_indices)
            other_ring_atoms = set(ring_indices) - set(unit_ring_indices)

            def add_substituents(atom_idx):
                atom = mol.GetAtomWithIdx(atom_idx)
                for neighbor in atom.GetNeighbors():
                    nb_idx = neighbor.GetIdx()
                    # Only add if neighbor isn't part of the core unit or the rest of the ring
                    if nb_idx not in final_atom_set and nb_idx not in other_ring_atoms:
                        final_atom_set.add(nb_idx)
                        add_substituents(nb_idx)

            for idx in unit_ring_indices:
                add_substituents(idx)

            # Create an editable molecule to add dummy atoms
            rw_mol = Chem.RWMol(mol)

            # Identify "Exit Bonds": Bonds from our unit to the rest of the ring
            exit_bonds = []
            for atom_idx in unit_ring_indices:
                atom = mol.GetAtomWithIdx(atom_idx)
                for neighbor in atom.GetNeighbors():
                    nb_idx = neighbor.GetIdx()
                    # If neighbor is in the ring but NOT in our current repeating unit
                    if nb_idx in other_ring_atoms:
                        exit_bonds.append((atom_idx, nb_idx))

            # Add dummy atoms at the exit points
            for stay_idx, leave_idx in exit_bonds:
                dummy_idx = rw_mol.AddAtom(Chem.Atom(0))
                rw_mol.AddBond(
                    stay_idx,
                    dummy_idx,
                    mol.GetBondBetweenAtoms(stay_idx, leave_idx).GetBondType(),
                )

            # Remove all atoms NOT in our final_atom_set or dummy atoms
            atoms_to_keep = final_atom_set | {
                idx for idx in range(mol.GetNumAtoms(), rw_mol.GetNumAtoms())
            }
            # remove in reverse to avoid indexing issues
            for idx in reversed(range(rw_mol.GetNumAtoms())):
                if idx not in atoms_to_keep:
                    rw_mol.RemoveAtom(idx)

            return rw_mol, True

    return mol, False


def break_cyclic_bond(cyclic_mol):
    """
    Break cyclic bond and return as linear mol with dummy atoms

    :param cyclic_mol: mol object for cyclic mol
    :return: mol object for linear mol
    """
    # reset indexing
    cyclic_mol = Chem.MolFromSmiles(Chem.MolToSmiles(cyclic_mol))
    cyclic_mol = Chem.RWMol(cyclic_mol)

    # Identify bonds in backbone
    for b in cyclic_mol.GetBonds():
        if b.IsInRing():
            b.SetProp("backbone", "True")

    # Find bonds to break
    bond_to_remove = []
    for b in cyclic_mol.GetBonds():
        if b.HasProp("backbone") and not b.HasProp("ring"):
            bond_to_remove.append((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
            break

    # Break cyclic bond and add dummy atoms back
    for begin_idx, end_idx in bond_to_remove:
        cyclic_mol.RemoveBond(begin_idx, end_idx)

        # Add first dummy atom (*) and bond it to begin_idx
        new_idx_1 = cyclic_mol.AddAtom(Chem.Atom(0))
        cyclic_mol.AddBond(begin_idx, new_idx_1, Chem.BondType.SINGLE)

        # Add second dummy atom (*) and bond it to end_idx
        new_idx_2 = cyclic_mol.AddAtom(Chem.Atom(0))
        cyclic_mol.AddBond(end_idx, new_idx_2, Chem.BondType.SINGLE)

    mol = cyclic_mol.GetMol()
    Chem.SanitizeMol(mol)
    return mol


def canonicalise(smiles):
    """
    Canonicalise a polymer smiles string and reduce multiples of a repeating unit

    :param smiles: polymer smiles string
    :return: canonical polymer smiles string
    """
    mol = Chem.MolFromSmiles(smiles)

    # Identify the dummy atoms (*) and their neighbors (anchors)
    dummy_indices = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "*"]
    anchors = [
        mol.GetAtomWithIdx(idx).GetNeighbors()[0].GetIdx() for idx in dummy_indices
    ]

    # check for only one atom in backbone e.g *C*
    if anchors[0] == anchors[1]:
        return Chem.MolToSmiles(mol)

    # check for bond between anchor atoms, e.g *CC* - these cant be cyclised
    if mol.GetBondBetweenAtoms(anchors[0], anchors[1]):
        return identify_symmetry(mol, anchors)

    # Identify bonds in rings
    for b in mol.GetBonds():
        if b.IsInRing():
            b.SetProp("ring", "True")

    # check if all bonds are in rings - these shouldnt be cyclised
    if all_ring_bonds(mol):
        return Chem.MolToSmiles(mol)

    # Create the cyclic version
    cyclic_mol = Chem.RWMol(mol)
    cyclic_mol.AddBond(anchors[0], anchors[1], Chem.rdchem.BondType.SINGLE)

    # Remove dummies in reverse order to keep indices stable
    for idx in sorted(dummy_indices, reverse=True):
        cyclic_mol.RemoveAtom(idx)

    # reduce multiplication
    cyclic_mol, repeating = get_repeating_substructure(cyclic_mol.GetMol())
    if repeating:
        return Chem.MolToSmiles(cyclic_mol)

    # return to linear form
    mol = break_cyclic_bond(cyclic_mol)
    return Chem.MolToSmiles(mol)


def find_canonical_repeats(smiles: str) -> list[str] or str:
    """
    Find repeat unit, and canonicalise a polymer smiles string.
    REMOVES END GROUPS
    Args:
        smiles - the unedited smiles from ketcher output. e.g. CC{-}C{+n}C
    Returns:
        smiles_list - list of the canonicalised repeat units. e.g ['*C*']
    """
    if "[*]" in smiles:  # block dummy atoms like R groups
        return "dummy"

    repeat_units = find_polymer_repeat_units(smiles)

    smiles_list = []
    for unit in repeat_units:
        unit = unit.replace("/", "").replace("\\", "")  # ignore stereochemistry
        smiles_list.append(canonicalise(unit))

    return smiles_list
