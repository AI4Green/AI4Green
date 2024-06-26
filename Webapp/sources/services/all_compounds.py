from typing import List, Optional, Tuple, Union

from rdkit import Chem
from sources import models, services
from sources.extensions import db

"""
Contains functions for when it is unknown if a compound is from the Compound table and from PubChem or a Novel compound
which is from the NovelCompound table and associated with a specific workgroup
"""


def get_smiles_list(primary_key_ls: List[Union[Tuple[str, int], int]], person: models.Person = None) -> List[str]:
    """
    Gets SMILES of compounds that could be type Compound (pk is int) or type NovelCompound (pk is Tuple[str, int])

    Args:
        primary_key_ls - the list of primary keys of the compounds
        person - the person we are checking for permission if accessing novel compounds

    Returns:
        The list of smiles for the compounds - or None if that item has no SMILES.
    """
    smiles_ls = []
    for primary_key in primary_key_ls:
        if validate_primary_key(primary_key):
            primary_key = primary_key_resolver(primary_key)
            smiles_ls.append(smiles_from_primary_key(primary_key, person=person))
    return smiles_ls


def smiles_from_primary_key(primary_key: Union[int, Tuple, str], person: models.Person = None) -> str:
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
