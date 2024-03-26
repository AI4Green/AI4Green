from typing import List, Tuple, Union

from sources import services

"""
Contains functions for when it is unknown if a compound is from the Compound table and from PubChem or a Novel compound
which is from the NovelCompound table and associated with a specific workgroup
"""


def get_smiles_list(primary_key_ls: List[Union[Tuple[str, int], int]]) -> List[str]:
    """
    Gets SMILES of compounds that could be type Compound (pk is int) or type NovelCompound (pk is Tuple[str, int])

    Args:
        primary_key_ls - the list of primary

    Returns:
        The list of smiles for the compounds - or None if that item has no SMILES.
    """
    smiles_ls = []
    for primary_key in primary_key_ls:
        if validate_primary_key(primary_key):
            primary_key = primary_key_resolver(primary_key)
            smiles_ls.append(smiles_from_primary_key(primary_key))
    return smiles_ls


def smiles_from_primary_key(primary_key: Union[int, Tuple, str]) -> str:
    if isinstance(primary_key, int):
        smiles = services.compound.get_smiles(primary_key)
    elif isinstance(primary_key, tuple):
        smiles = services.novel_compound.get_smiles(primary_key)
    return smiles


def primary_key_resolver(primary_key: Union[int, Tuple, str]) -> Union[int, Tuple]:
    """
    In the compound table int is the primary key and in the novel compound table tuple is the primary key
    The frontend may send these back as a string. This function converts it to either an int or a tuple
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


def validate_primary_key(primary_key):
    """Validates primary key is real value and not a default database entry value or None"""
    if primary_key:
        if primary_key.isdigit():
            return True
        if len(primary_key) > 1:
            return True
    return False
