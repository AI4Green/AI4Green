import os
import pickle
import tempfile

import chython
import pytest
import pytest_mock
from flask import Flask
from sources import services

IN_GITHUB_ACTIONS = os.getenv("GITHUB_ACTIONS") == "true"


@pytest.mark.skipif(IN_GITHUB_ACTIONS, reason="Local temp files don't work in CI")
def test_export_reaction_as_rdf(app: Flask, mocker: pytest_mock.MockerFixture):
    """
    Test saving an exported reaction data file locally and reloading it and checking for data changes

    This test ensures that the services.data_export.ReactionDataFile.save_as_rdf function correctly saves the .RDF
    It then tests we can reload it and the contents are the same as the original.

    Args:
        app: Flask test client.
        mocker: Mock

    """
    with open(
        os.path.join(
            os.path.dirname(os.getcwd()), "data", "reaction_database_object.pickle"
        ),
        "rb",
    ) as f:
        serialized_data = f.read()
    # mock the database responses
    make_mocks(mocker)

    reaction = pickle.loads(serialized_data)
    with app.app_context():
        # making a temporary directory to save the output file
        local_path = tempfile.gettempdir()
        file_path = os.path.join(local_path, "test_rdf")
        # take a reaction pickle and confirm the function exports correctly.
        test_rdf = services.data_export.ReactionDataFile(reaction, file_path)

        # save and load the file and confirm no data has been lost
        rdf = test_rdf.reaction_container
        test_rdf.save_as_rdf()

        with chython.files.RDFRead(file_path) as f:
            rdf_contents = next(f)
        test_rdf.literal_eval_metadata(rdf_contents)
        validate_rdf(rdf, rdf_contents)


def make_mocks(mocker: pytest_mock.MockerFixture):
    """Makes mocks and return values for the functions that would normally access the database during RDF creation."""
    mock_get_smiles_from_db = mocker.patch(
        "sources.services.all_compounds.smiles_from_primary_key"
    )
    mock_get_smiles_from_db.return_value = "CO"

    mock_read_username = mocker.patch(
        "sources.services.data_export.ReactionMetaData.read_creator_username"
    )
    mock_read_username.return_value = "Al Forgreen"

    mock_read_workbook = mocker.patch(
        "sources.services.data_export.ReactionMetaData.read_workbook"
    )
    mock_read_workbook.return_value = "Alchemical aspirations"

    mock_read_workgroup = mocker.patch(
        "sources.services.data_export.ReactionMetaData.read_workgroup"
    )
    mock_read_workgroup.return_value = "AI4Alchemy"

    mock_read_file_attachments = mocker.patch(
        "sources.services.data_export.ReactionMetaData.read_file_attachment_names"
    )
    mock_read_file_attachments.return_value = ["AI1-001-summary.pdf"]

    mock_read_addenda = mocker.patch(
        "sources.services.data_export.ReactionMetaData.read_addenda"
    )
    mock_read_addenda.return_value = None


def validate_rdf(
    reaction_container: chython.ReactionContainer,
    rdf_contents: chython.ReactionContainer,
):
    """Validates the read RDF file contains the same data as the original object in Python"""
    assert rdf_contents.meta == reaction_container.meta, "unequal  metas"

    for original_reactant, loaded_reactant in zip(
        rdf_contents.reactants, reaction_container.reactants
    ):
        if original_reactant and loaded_reactant:
            assert original_reactant.is_equal(loaded_reactant), "reactant inequality"

    for original_reagent, loaded_reagent in zip(
        rdf_contents.reagents, reaction_container.reagents
    ):
        if original_reagent and loaded_reagent:
            assert original_reagent.is_equal(loaded_reagent), "reagent inequality"

    for original_product, loaded_product in zip(
        rdf_contents.products, reaction_container.products
    ):
        if original_product and loaded_product:
            assert original_product.is_equal(loaded_product), "product inequality"
