"""
Chemical structure matching shared by the reaction search and the inventory.

A search is built once from the structure the user drew (`build_query`) and is
then compared with candidate SMILES strings (`matches` / `any_match`).

- Exact searches compare standard InChI strings, so different SMILES for the
  same molecule still match.
- Substructure searches use RDKit's `HasSubstructMatch`, so a candidate matches
  when it contains the drawn fragment.
- Polymers (SMILES with `*` end groups) are matched exactly by their
  canonicalised repeat unit, the same way the reaction search always has.
  Substructure search isn't supported for them.

This module has no database or Flask dependencies, so callers decide which
candidates the user is allowed to search.
"""

from dataclasses import dataclass
from enum import Enum
from functools import lru_cache
from typing import Iterable, Optional

from rdkit import Chem
from sources import services

# Longest structure accepted from the sketcher. This is far beyond any real
# molecule, and it stops a huge string being parsed.
MAX_QUERY_SMILES_LENGTH = 5000
# How many parsed candidate structures to keep between searches
CANDIDATE_CACHE_SIZE = 4096


class SearchMode(str, Enum):
    """How the drawn structure is compared with each candidate."""

    EXACT = "exact"
    SUBSTRUCTURE = "substructure"


class StructureSearchError(ValueError):
    """The drawn structure can't be searched. The message is safe to show to the user."""


@dataclass(frozen=True)
class StructureQuery:
    """
    A structure ready to compare with candidates.

    Attributes:
        mode: Exact or substructure matching.
        smiles: The drawn SMILES without any CXSMILES extension.
        is_polymer: Whether this is a polymer repeat unit (contains `*`).
        inchi: The InChI to compare with, for exact non-polymer searches.
        canonical_polymer_smiles: The canonical repeat unit, for polymer searches.
        mol: The RDKit molecule to look for, for substructure searches.
    """

    mode: SearchMode
    smiles: str
    is_polymer: bool
    inchi: Optional[str] = None
    canonical_polymer_smiles: Optional[str] = None
    mol: Optional[Chem.Mol] = None


def read_mode(value: object) -> SearchMode:
    """
    Reads a search mode sent by the front end.

    Args:
        value: "exact" or "substructure".

    Returns:
        The search mode.

    Raises:
        StructureSearchError: If the value isn't a known mode.
    """
    try:
        return SearchMode(value)
    except ValueError:
        raise StructureSearchError(
            "Choose either an exact or a substructure search."
        ) from None


def strip_cxsmiles(smiles: str) -> str:
    """
    Removes the CXSMILES extension (everything after " |") that Ketcher adds,
    e.g. for polymer brackets, leaving plain SMILES.

    Args:
        smiles: A SMILES or CXSMILES string.

    Returns:
        The plain SMILES.
    """
    return smiles.split(" |")[0].strip()


def build_query(smiles: object, mode: SearchMode) -> StructureQuery:
    """
    Checks the drawn structure and prepares it for matching.

    Args:
        smiles: The SMILES exported from the sketcher.
        mode: Exact or substructure.

    Returns:
        The prepared query.

    Raises:
        StructureSearchError: If there's no structure, it's a reaction, it
            can't be parsed, or it's a polymer substructure search.
    """
    if not isinstance(smiles, str) or not smiles.strip():
        raise StructureSearchError("Draw a structure to search for.")
    if len(smiles) > MAX_QUERY_SMILES_LENGTH:
        raise StructureSearchError("The structure is too large to search for.")

    plain_smiles = strip_cxsmiles(smiles)
    if ">" in plain_smiles:
        raise StructureSearchError(
            "Draw a single structure to search for, not a reaction."
        )

    if "*" in plain_smiles:
        return _build_polymer_query(plain_smiles, mode)

    mol = Chem.MolFromSmiles(plain_smiles)
    if mol is None:
        raise StructureSearchError("Invalid structure, please try again!")

    if mode is SearchMode.EXACT:
        return StructureQuery(
            mode=mode,
            smiles=plain_smiles,
            is_polymer=False,
            inchi=Chem.MolToInchi(mol),
        )
    return StructureQuery(mode=mode, smiles=plain_smiles, is_polymer=False, mol=mol)


def _build_polymer_query(plain_smiles: str, mode: SearchMode) -> StructureQuery:
    """
    Prepares a polymer repeat unit for an exact search.

    Args:
        plain_smiles: The repeat unit's SMILES, with `*` end groups.
        mode: Exact or substructure.

    Returns:
        The prepared query.

    Raises:
        StructureSearchError: For substructure searches, polymers with known end
            groups, or SMILES that can't be parsed.
    """
    if mode is SearchMode.SUBSTRUCTURE:
        raise StructureSearchError(
            "Substructure search is not yet supported for polymers."
        )
    if plain_smiles.count("*") != 2:
        raise StructureSearchError(
            "Searching for polymers with known end groups is not yet supported."
        )
    if Chem.MolFromSmiles(plain_smiles) is None:
        raise StructureSearchError("Invalid structure, please try again!")

    canonical = services.polymer_novel_compound.canonicalise(plain_smiles)
    return StructureQuery(
        mode=mode,
        smiles=plain_smiles,
        is_polymer=True,
        canonical_polymer_smiles=canonical,
    )


@lru_cache(maxsize=CANDIDATE_CACHE_SIZE)
def _mol_from_smiles(smiles: str) -> Optional[Chem.Mol]:
    """
    Parses a candidate SMILES, caching the result so compounds that appear in
    many reactions are only parsed once.

    Args:
        smiles: Plain candidate SMILES.

    Returns:
        The molecule, or None if the SMILES can't be parsed.
    """
    return Chem.MolFromSmiles(smiles)


@lru_cache(maxsize=CANDIDATE_CACHE_SIZE)
def _inchi_from_smiles(smiles: str) -> Optional[str]:
    """
    Gets a candidate's InChI, caching the result.

    Args:
        smiles: Plain candidate SMILES.

    Returns:
        The InChI, or None if the SMILES can't be parsed.
    """
    mol = _mol_from_smiles(smiles)
    return None if mol is None else Chem.MolToInchi(mol)


def matches(query: StructureQuery, candidate_smiles: Optional[str]) -> bool:
    """
    Whether a candidate structure matches the query.

    Args:
        query: The prepared query.
        candidate_smiles: The candidate's SMILES. Empty or invalid SMILES never match.

    Returns:
        True if the candidate is the queried structure (exact) or contains it
        (substructure).
    """
    if not candidate_smiles or not isinstance(candidate_smiles, str):
        return False
    plain_smiles = strip_cxsmiles(candidate_smiles)
    if not plain_smiles:
        return False

    if query.is_polymer:
        return (
            services.polymer_novel_compound.clean_polymer_smiles(plain_smiles)
            == query.canonical_polymer_smiles
        )
    if query.mode is SearchMode.EXACT:
        return _inchi_from_smiles(plain_smiles) == query.inchi

    candidate_mol = _mol_from_smiles(plain_smiles)
    return candidate_mol is not None and candidate_mol.HasSubstructMatch(query.mol)


def any_match(query: StructureQuery, candidates: Iterable[Optional[str]]) -> bool:
    """
    Whether any of the candidate structures matches the query.

    Args:
        query: The prepared query.
        candidates: Candidate SMILES strings.

    Returns:
        True if at least one candidate matches.
    """
    return any(matches(query, smiles) for smiles in candidates)
