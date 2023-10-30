from typing import List, Tuple, Union

from sources import models
from sources.extensions import db
from sources.services import queries

"""
Contains functions for when it is unknown if a compound is from the Compound table and from PubChem or a Novel compound
which is from the NovelCompound table and associated with a specific workgroup
"""


def get_smiles_from_primary_keys(
    primary_key_ls: List[Union[Tuple[str, int], int]]
) -> List[str]:
    """
    Gets SMILES of compounds that could be type Compound (pk is int) or type NovelCompound (pk is Tuple[str, int])

    Args:
        primary_key_ls - the list of primary

    Returns:
        The list of smiles for the compounds - or None if that item has no SMILES.
    """
    smiles_ls = []
    for primary_key in primary_key_ls:
        if isinstance(primary_key, int):
            smiles = queries.compound.get_smiles_from_primary_key(primary_key)
        else:
            smiles = queries.novel_compound.get_smiles_from_primary_key(primary_key)
        smiles_ls.append(smiles)
    return smiles_ls
