from typing import List, Optional, Tuple, Union

from rdkit import Chem
from rdkit.Chem import Descriptors
from sources import models, services
from sources.extensions import db
from sqlalchemy import func

"""
Contains functions for when it is unknown if a compound is from the Compound table and from PubChem or a Novel compound
which is from the NovelCompound table and associated with a specific workgroup
"""


def get_smiles_list(
    primary_key_ls: List[Union[Tuple[str, int], int, str]],
    polymer_indices: List[str] = None,
    number_of_reactants: int = 0,
    person: models.Person = None,
) -> List[str]:
    """
    Gets SMILES of compounds that could be type Compound (pk is int) or type NovelCompound (pk is Tuple[str, int])

    Args:
        primary_key_ls - the list of primary keys of the compounds
        polymer_indices - list of indices where polymers are in reaction
        number_of_reactants - used to update product indexing to check for polymers, can be empty when checking reactant smiles
        person - the person we are checking for permission if accessing novel compounds

    Returns:
        The list of smiles for the compounds - or None if that item has no SMILES.
    """
    smiles_ls = []

    if polymer_indices is None:
        polymer_indices = []

    for idx, primary_key in enumerate(primary_key_ls):
        if validate_primary_key(primary_key):
            primary_key = primary_key_resolver(primary_key)

            if (idx + 1 + number_of_reactants) in polymer_indices:
                smiles_ls.append(
                    services.polymer_novel_compound.get_smiles(  # get smiles from polymer db
                        primary_key, person=person
                    )
                )
            else:
                smiles_ls.append(smiles_from_primary_key(primary_key, person=person))
    return smiles_ls


def smiles_from_primary_key(
    primary_key: Union[int, Tuple, str], person: models.Person = None
) -> str:
    """
    Gets the SMILES for a compound or a novel compound from the primary key

    Args:
        primary_key - integer for compound or a Tuple for novel compound
        person - the person we are checking for permission to access the novel compounds
    Returns:
        SMILEs string for the compound or novel compound

    """
    if isinstance(primary_key, int):
        smiles = services.compound.get_smiles(primary_key)
    elif isinstance(primary_key, tuple):
        smiles = services.novel_compound.get_smiles(primary_key, person=person)
    return smiles


def primary_key_resolver(primary_key: Union[int, Tuple, str]) -> Union[int, Tuple]:
    """
    In the compound table int is the primary key and in the novel compound table tuple is the primary key
    The frontend may send these back as a string. This function converts it to either an int or a tuple

    Args:
        primary_key - either an int, tuple, or the string equivalent of one of these
    Returns:
        The primary key as an integer for compound or tuple for novel compound

    """
    if primary_key.isdigit():
        primary_key = int(primary_key)
    elif isinstance(primary_key, tuple):
        primary_key = primary_key
    else:
        primary_key = services.novel_compound.reform_novel_compound_primary_key(
            primary_key
        )
    return primary_key


def validate_primary_key(primary_key: Union[int, Tuple, str]) -> bool:
    """
    Validates primary key is real value and not a default database value or None

    Args:
        primary_key - either an int, tuple, or the string equivalent of one of these

    Returns:
        True if the primary key is valid
    """
    if primary_key:
        if primary_key.isdigit():
            return True
        if len(primary_key) > 1:
            return True
    return False


def cas_from_smiles(smiles: str) -> Optional[str]:
    """
    Look in compound and novel compound database. Return CAS if present in either else None

    Args:
        smiles - the SMILES string of the compound we want the CAS number for

    Returns:
        The CAS number of the compound or None if no cas number is provided


    """
    cas = None
    mol = Chem.MolFromSmiles(smiles)
    inchi = Chem.MolToInchi(mol)
    compound = (
        db.session.query(models.Compound).filter(models.Compound.inchi == inchi).first()
    )
    if compound:
        cas = compound.cas
    else:
        novel_compound = (
            db.session.query(models.Compound)
            .filter(models.NovelCompound.inchi == inchi)
            .first()
        )
        if novel_compound:
            cas = novel_compound.cas
    return cas if cas else None


def name_from_smiles(smiles: str) -> Optional[str]:
    """

    Look in compound and novel compound database. Return Name if present in either else None

    Args:
        smiles - the SMILES string of the compound we want the CAS number for

    Returns:
        The name of the compound or None if no name is provided


    """
    name = None
    compound = services.compound.from_smiles(smiles)
    if compound:
        name = compound.name
    else:
        inchi = smiles_to_inchi(smiles)
        novel_compound = (
            db.session.query(models.Compound)
            .filter(models.NovelCompound.inchi == inchi)
            .first()
        )
        if novel_compound:
            name = novel_compound.name
    return name if name else None


def smiles_to_inchi(smiles: str) -> str:
    """
    Get the corresponding InChi from a SMILES string. Returns None if the smiles string is invalid

    Args:
       smiles: The compound's smiles string

    Returns:
        The compound's InChI string
    """
    mol = Chem.MolFromSmiles(smiles)
    return None if mol is None else Chem.MolToInchi(mol)


def mol_weight_from_smiles(smiles: Union[str, List]) -> Union[float, List]:
    """
    Uses RDKit to calculate the molecular weight for a compound from its SMILES string

    Args:
        smiles - the SMILES of the compound of interest

    Returns:
        The molecular weight of the compound.
    """
    # MolWt accounts for the average across isotopes but ExactMolWt only takes the most abundant isotope.
    if isinstance(smiles, list):
        mol_wts = [round(Descriptors.MolWt(Chem.MolFromSmiles(s)), 2) for s in smiles]
    else:
        mol_wts = round(Descriptors.MolWt(Chem.MolFromSmiles(smiles)), 2)

    return mol_wts


def from_cas(
    cas: str, workbook: models.WorkBook = None
) -> Union[models.Compound, models.NovelCompound]:
    """
    Returns a compound or novel compound from a cas string

    Args:
        cas - cas number of a compound
        workbook - if looking for a novel compound, the workbook we look in

    Returns:
        models.Compound or models.NovelCompound: The compound object retrieved from the database.
                                                 Returns None if no matching compound is found.
    """
    compound = services.compound.from_cas(cas)
    if workbook and not compound:
        compound = services.novel_compound.from_cas_and_workbook(cas, workbook)
    return compound


def from_name(
    name: str, workbook: models.WorkBook = None
) -> Union[models.Compound, models.NovelCompound]:
    """
    Returns a compound or novel compound from a name

    Args:
        name - name of a compound
        workbook - if looking for a novel compound, the workbook we look in

    Returns:
        models.Compound or models.NovelCompound: The compound object retrieved from the database.
                                                 Returns None if no matching compound is found.
    """
    compound = services.compound.from_name(name)
    if workbook and not compound:
        compound = services.novel_compound.from_name_and_workbook(name, workbook)
    return compound


def from_inchi(
    inchi: str, workbook: models.WorkBook = None
) -> Union[models.Compound, models.NovelCompound]:
    """
    Returns a compound or novel compound from an inchi

    Args:
        inchi - inchi of a compound
        workbook - if looking for a novel compound, the workbook we look in

    Returns:
        models.Compound or models.NovelCompound: The compound object retrieved from the database.
                                                 Returns None if no matching compound is found.
    """
    compound = services.compound.from_inchi(inchi)
    if workbook and not compound:
        compound = services.novel_compound.from_inchi_and_workbook(inchi, workbook)
    return compound


def populate_reagent_dropdown(
    reagent_substring: str, workbook: models.WorkBook = None
) -> List[str]:
    """
    Makes the dropdown for the reagent input field in the reaction constructor.
    When a user first clicks the reagent input, the substring will be an empty string
    and only novel compounds/recent reagents will populate the list.
    Once the user starts typing the list will filter by the substring for novel compounds/recent reagents first
    and the rest of the list up to 100 is completed by compounds from the Compound Table.

    Args:
        reagent_substring - used to filter by substring. An empty substring is ignored
        workbook - the active workbook
    Returns:
        A list of names of novel compounds/compounds.

    """
    remaining_spaces = 100
    reagent_names = []
    # if the user is making a reaction in a workbook, get its novel compounds and recent reagents
    if workbook:
        novel_compound_list = services.novel_compound.all_from_workbook(workbook)
        full_reagent_list = novel_compound_list + workbook.recent_compounds
        reagent_list = services.utils.remove_duplicates_keep_first(full_reagent_list)

        # Get the first 100 which match the substring
        reagent_names = [
            x.name for x in reagent_list if reagent_substring.lower() in x.name.lower()
        ][0:100]
        # Calculate remaining spaces needed to fill up to 100
        remaining_spaces -= len(reagent_names)

    # Query for additional reagent from the Compounds Table to fill up the list to 100
    if remaining_spaces > 0 and reagent_substring:
        additional_reagents = (
            db.session.query(models.Compound)
            .filter(
                func.lower(models.Compound.name).startswith(reagent_substring.lower())
            )
            .order_by(models.Compound.name.asc())
            .limit(remaining_spaces)
            .all()
        )
        reagent_names.extend([x.name for x in additional_reagents])
    return reagent_names


def check_valid_smiles(smiles: str or list):
    """
    Checks if a smiles string or list of smiles strings are valid.

    Args:
        smiles: The compound's smiles string

    Returns:
        True if compound is valid, False otherwise.
    """
    if isinstance(smiles, list):
        for s in smiles:
            mol = Chem.MolFromSmiles(s)
            if mol is None:
                return False
    else:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False

    return True
