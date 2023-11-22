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
        if isinstance(primary_key, int):
            smiles = services.compound.get_smiles(primary_key)
        else:
            smiles = services.novel_compound.get_smiles(primary_key)
        smiles_ls.append(smiles)
    return smiles_ls
