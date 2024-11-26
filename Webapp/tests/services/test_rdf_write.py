import ast
import os
import pickle
import tempfile
from pathlib import Path
from typing import Dict, List

import pytest
import pytest_mock
from flask import Flask
from flask.testing import FlaskClient
from rdkit import Chem
from rdkit.Chem import AllChem
from sources import db, services

IN_GITHUB_ACTIONS = os.getenv("GITHUB_ACTIONS") == "true"


@pytest.mark.skipif(IN_GITHUB_ACTIONS, reason="Local temp files don't work in CI")
def test_export_reaction_as_rdf(
    client: FlaskClient, app: Flask, mocker: pytest_mock.MockerFixture
):
    """
    Test saving an exported reaction data file locally and reloading it and checking for data changes

    This test ensures that the services.data_export.ReactionDataFile.save_as_rdf function correctly saves the .RDF
    It then tests we can reload it and the contents are the same as the original.

    Args:
        app: Flask test client.
        mocker: Mock

    """
    with open(
        Path(__file__).resolve().parent.parent
        / "data"
        / "reaction_database_object_v1_5.pickle",
        "rb",
    ) as f:
        serialized_data = f.read()
    # mock the database responses
    make_mocks(mocker)

    reaction = pickle.loads(serialized_data)
    with app.app_context():
        session = db.session()
        session.add(reaction)
        # making a temporary directory to save the output file
        local_path = tempfile.gettempdir()
        file_path = os.path.join(local_path, "test_rdf")
        # take a reaction pickle and confirm the function exports correctly.
        test_rdf = services.data_export.serialisation.ReactionDataFileExport(
            reaction, file_path, "containertest"
        )
        rdf_contents = test_rdf._make_rdf_contents()

        # save and load the file and confirm no data has been lost
        saved_file_path = os.path.join(local_path, "test_rdf.rdf")
        with open(saved_file_path, "w") as f:
            f.write(rdf_contents)

        # read rdf
        previous_line = ""
        metadata = {}
        key = ""
        in_rxn_block = False
        rxn_block = ""

        with open(saved_file_path, "r") as f:
            for line in f:
                # get the role indices
                if previous_line == "RDKit":
                    role_indices = line.strip().split("  ")
                # get the rxnblock data
                if line.strip() == "$RXN":
                    in_rxn_block = True

                if line.startswith("$DYPTE"):
                    in_rxn_block = False

                if in_rxn_block is True:
                    rxn_block += line
                # get the metadata
                if line.startswith("$DTYPE"):
                    key = line.strip().split("$DTYPE ")[-1]
                if previous_line.startswith("$DTYPE "):
                    metadata.update({key: line.strip().split("$DATUM ")[-1]})
                previous_line = line.strip() if line.strip() else previous_line

        role_indices = [int(x) for x in role_indices]
        literal_eval_metadata(metadata)
        validate_rdf(test_rdf, role_indices, rxn_block, metadata)


def literal_eval_metadata(metadata: Dict):
    """Read the rdf values literally to reintroduce their types"""
    non_string_types_inside_metadata = [
        "reagents",
        "reactants",
        "solvents",
        "products",
        "standard_protocols_used",
        "file_attachment_names",
        "addenda",
        "solvent_sustainability",
        "sustainability_data",
        "polymer_data",
    ]

    for key in metadata.keys():
        print(key)
        print(metadata[key])
        # if the string is equal to none we reload
        if (
            metadata[key] == "None"
            or key in non_string_types_inside_metadata
            and isinstance(metadata[key], str)
        ):
            metadata.update({key: ast.literal_eval(metadata[key])})


def make_mocks(mocker: pytest_mock.MockerFixture):
    """Makes mocks and return values for the functions that would normally access the database during RDF creation."""
    mock_get_smiles_from_db = mocker.patch(
        "sources.services.all_compounds.smiles_from_primary_key"
    )
    mock_get_smiles_from_db.return_value = "CO"

    mock_read_username = mocker.patch(
        "sources.services.data_export.metadata.ReactionMetaData._read_creator_username"
    )
    mock_read_username.return_value = "Al Forgreen"

    mock_read_workbook = mocker.patch(
        "sources.services.data_export.metadata.ReactionMetaData._read_workbook"
    )
    mock_read_workbook.return_value = "Alchemical aspirations"

    mock_read_workgroup = mocker.patch(
        "sources.services.data_export.metadata.ReactionMetaData._read_workgroup"
    )
    mock_read_workgroup.return_value = "AI4Alchemy"

    mock_read_file_attachments = mocker.patch(
        "sources.services.data_export.metadata.ReactionMetaData._read_file_attachment_names"
    )
    mock_read_file_attachments.return_value = ["AI1-001-summary.pdf"]

    mock_read_addenda = mocker.patch(
        "sources.services.data_export.metadata.ReactionMetaData._read_addenda"
    )
    mock_read_addenda.return_value = None


def validate_rdf(
    original_rdf: services.data_export.serialisation.ReactionDataFileExport,
    role_indices: List[int],
    rxn_block: str,
    metadata: Dict,
):
    """
    Validates the read RDF file contains the expected data
    Args:
        original_rdf - the rdf that was made during the export
        role_indices - the number of each reactant, product and agents (in that order)
        rxn_block - the rxn block read from the rdf file
        metadata - the metadata read from the rdf file and literally evaluated
    """
    assert role_indices == [2, 1, 2], "unequal role indices"

    rxn = AllChem.ReactionFromRxnBlock(rxn_block)
    rxn_smiles = AllChem.ReactionToSmiles(rxn, canonical=True)
    reactants = [Chem.CanonSmiles(smi) for smi in rxn_smiles.split(">")[0].split(".")]
    agents = [Chem.CanonSmiles(smi) for smi in rxn_smiles.split(">")[1].split(".")]
    products = [Chem.CanonSmiles(smi) for smi in rxn_smiles.split(">")[-1].split(".")]

    # combine reagent and solvent SMILES to get agents
    original_agents = (
        original_rdf.db_reaction.reagents
        + services.all_compounds.get_smiles_list(original_rdf.db_reaction.solvent)
    )
    original_agents = [x for x in original_agents if x]

    assert original_rdf.metadata == metadata

    assert sorted(reactants) == sorted(
        original_rdf.db_reaction.reactants
    ), "difference in reactants"

    assert sorted(agents) == sorted(original_agents), "difference in reagents"

    assert sorted(products) == sorted(
        original_rdf.db_reaction.products
    ), "difference in products"
