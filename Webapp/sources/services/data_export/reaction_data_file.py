import ast
import io
import json

import chython.files
from flask import abort
from flask_login import current_user
from sources import models, services

from . import utils


class ReactionDataFile:
    non_string_types_inside_meta = [
        "reagents",
        "reactants",
        "solvents",
        "products",
        "standard_protocols_used",
        "file_attachment_names",
        "addenda",
        "solvent_sustainability",
        "sustainability_data",
    ]

    def __init__(
        self, db_reaction: models.Reaction, filename: str, container_name: str
    ):
        """
        Filename should not include extension.
        """
        self.db_reaction = db_reaction
        self.reaction_container = self._make_reaction_container()
        self.metadata = services.data_export.metadata.ReactionMetaData(
            db_reaction
        ).get_dict()
        self.container_name = container_name
        self.filename = filename
        self.memory_file = None
        self.file_hash = None
        if self.reaction_container:
            self.reaction_container.meta.update(self.metadata)

    def save(self):
        """Calls the appropriate method to save the data. Can't use .RDF without a reaction_container object"""
        if self.reaction_container:
            self.filename += ".rdf"
            self._save_as_rdf()
        else:
            self.filename += ".json"
            self._save_as_json()
        self._save_blob()

    def _make_reaction_container(self):
        """
        Makes a reaction container by using the SMILES obtained from the db reaction object.
        Full reaction SMILES including reactants, agents (aka reagents and solvents), and products
        """
        reaction_smiles = self._make_reaction_smiles()
        # return None if we cannot make a reaction_container from smiles - likely to be due to no smiles present.
        try:
            return chython.files.smiles(reaction_smiles)
        except ValueError:
            return None

    def _make_reaction_smiles(self) -> str:
        """
        Makes a reaction smiles in format reactant1.reactants2>reagents.solvents>products from a db reaction object

        Returns:
            the reaction smiles string.
        """
        reactant_smiles = utils.remove_default_data(self.db_reaction.reactants)
        reagent_smiles = utils.remove_default_data(self.db_reaction.reagents)
        solvent_smiles = services.all_compounds.get_smiles_list(
            self.db_reaction.solvent
        )
        product_smiles = utils.remove_default_data(self.db_reaction.products)
        reaction_smiles = (
            ".".join(reactant_smiles)
            + ">"
            + ".".join(reagent_smiles)
            + ".".join(solvent_smiles)
            + ">"
            + ".".join(product_smiles)
        )
        return reaction_smiles

    def _save_blob(self):
        # connect to azure
        blob_service_client = services.file_attachments.connect_to_azure_blob_service()
        # make more get container
        container_name = "data_export_temp"
        # container_client = (
        #     services.file_attachments.create_or_get_existing_container_client(
        #         blob_service_client, container_name
        #     )
        # )
        blob_client = blob_service_client.get_blob_client(
            container=container_name, blob=self.filename
        )

        print("uploading blob")
        # self.memory_file.stream.seek(0)
        # self.file_hash = services.file_attachments.sha256_from_file_contents(
        #     self.memory_file.stream.read()
        # )
        # if not self.file_hash:
        #     print("Failed to generate hash")
        #     abort(500)
        # Upload the blob data - default blob type is BlockBlob
        self.memory_file.stream.seek(0)
        upload = io.BytesIO(self.memory_file.stream.read())
        blob_client.upload_blob(upload, blob_type="BlockBlob")
        # confirm upload
        if not blob_client.exists():
            print(f"blob {self.filename} upload failed")
            abort(401)

    def _save_as_rdf(self):
        """Saves as an RDF"""
        # Create an in-memory file-like object
        self.memory_file = io.BytesIO()
        with chython.files.RDFWrite(self.filename) as f:
            f.write(self.reaction_container)

    def _save_as_json(self):
        with open(self.filename, "w") as f:
            json.dump(self.metadata, f)

        with open(self.filename, "r") as f:
            meta_reload = json.load(f)

        assert meta_reload == self.metadata, "change in data during file read/write"
