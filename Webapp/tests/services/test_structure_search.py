"""
Unit tests for services.structure_search, the RDKit matching shared by the
reaction search and the inventory. These need no database.
"""

import pytest
from sources.services.structure_search import (
    SearchMode,
    StructureSearchError,
    any_match,
    build_query,
    matches,
    read_mode,
)

BENZENE = "c1ccccc1"
TOLUENE = "Cc1ccccc1"
ACETIC_ACID = "CC(=O)O"
CARBOXYLIC_ACID = "C(=O)O"


def test_read_mode_accepts_known_modes():
    """The front end's mode strings are read into SearchMode values."""
    assert read_mode("exact") is SearchMode.EXACT
    assert read_mode("substructure") is SearchMode.SUBSTRUCTURE


@pytest.mark.parametrize("mode", [None, "", "similarity", 1])
def test_read_mode_rejects_unknown_modes(mode):
    """Anything else is a user-facing error."""
    with pytest.raises(StructureSearchError):
        read_mode(mode)


def test_exact_match_ignores_how_the_smiles_is_written():
    """Exact searches compare InChI, so equivalent SMILES match."""
    query = build_query(TOLUENE, SearchMode.EXACT)
    assert matches(query, "c1ccccc1C")
    assert matches(query, "CC1=CC=CC=C1")
    assert not matches(query, BENZENE)


def test_substructure_matches_candidates_containing_the_query():
    """Benzene is in toluene, but toluene isn't in benzene."""
    benzene = build_query(BENZENE, SearchMode.SUBSTRUCTURE)
    assert matches(benzene, TOLUENE)
    assert matches(benzene, BENZENE)
    assert not matches(benzene, "C1CCCCC1")  # cyclohexane isn't aromatic

    toluene = build_query(TOLUENE, SearchMode.SUBSTRUCTURE)
    assert not matches(toluene, BENZENE)


def test_substructure_finds_a_functional_group():
    """A carboxylic acid group is found in acetic acid but not in ethanol."""
    query = build_query(CARBOXYLIC_ACID, SearchMode.SUBSTRUCTURE)
    assert matches(query, ACETIC_ACID)
    assert not matches(query, "CCO")


def test_cxsmiles_extensions_are_ignored():
    """Ketcher may add a CXSMILES extension after " |", which isn't part of the structure."""
    cxsmiles = "Cc1ccccc1 |c:1,3,5|"
    query = build_query(cxsmiles, SearchMode.EXACT)
    assert query.smiles == TOLUENE
    assert matches(query, cxsmiles)


@pytest.mark.parametrize("candidate", [None, "", "not a smiles", "C1CC"])
def test_empty_or_invalid_candidates_never_match(candidate):
    """Bad stored SMILES are skipped rather than raising."""
    assert not matches(build_query("C", SearchMode.EXACT), candidate)
    assert not matches(build_query("C", SearchMode.SUBSTRUCTURE), candidate)


def test_any_match():
    """any_match is true if one candidate matches."""
    query = build_query(BENZENE, SearchMode.SUBSTRUCTURE)
    assert any_match(query, ["CCO", TOLUENE])
    assert not any_match(query, ["CCO", None])
    assert not any_match(query, [])


@pytest.mark.parametrize(
    "smiles, message",
    [
        (None, "Draw a structure"),
        ("   ", "Draw a structure"),
        ("thisIsDefinitelyNotASmilesString", "Invalid structure"),
        ("CC>>CO", "not a reaction"),
        ("C" * 5001, "too large"),
    ],
)
def test_invalid_queries_are_rejected(smiles, message):
    """Queries that can't be searched raise an error with a user-facing message."""
    with pytest.raises(StructureSearchError, match=message):
        build_query(smiles, SearchMode.SUBSTRUCTURE)


def test_polymer_exact_search():
    """Polymer repeat units are matched by their canonical form."""
    query = build_query("*C*", SearchMode.EXACT)
    assert query.is_polymer
    assert matches(query, "*C*")
    assert not matches(query, "*CC(C)*")


def test_polymer_substructure_search_is_not_supported():
    """Substructure search isn't defined for polymer repeat units."""
    with pytest.raises(StructureSearchError, match="not yet supported for polymers"):
        build_query("*C*", SearchMode.SUBSTRUCTURE)


def test_polymers_with_known_end_groups_are_not_supported():
    """Only repeat units with exactly two attachment points can be searched."""
    with pytest.raises(StructureSearchError, match="known end groups"):
        build_query("*C(*)*", SearchMode.EXACT)
