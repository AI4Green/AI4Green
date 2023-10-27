from typing import List, Tuple, Union

from sources import models
from sources.extensions import db
from sources.services import queries


def get_ambiguous_compound_smiles_from_primary_keys(
    primary_key_ls: List[Union[Tuple[str, int]], int]
) -> List[str]:
    """
    Gets SMILES of compounds that could be type Compound (pk is int) or type NovelCompound (pk is Tuple[str, int])

    Returns:
        The list of smiles for the compounds - or None if that item has no SMILES.
    """
    smiles_ls = []
    for primary_key in primary_key_ls:
        if isinstance(primary_key, int):
            smiles = queries.compound.get_smiles_from_primary_key()
        else:
            smiles = queries.novel_compound.get_smiles_from_primary_key()
        smiles_ls.append(smiles)
