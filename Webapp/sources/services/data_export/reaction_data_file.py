import datetime
import io
import json
from typing import Optional

import azure.core.exceptions
import pytz
from flask import abort
from rdkit.Chem import AllChem
from sources import models, services

from . import utils


class ReactionDataFile:
    # currently only used in testing but in future will want to read files in
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
    ]

    def __init__(
        self, db_reaction: models.Reaction, filename: str, container_name: str
    ):
        """
        Filename should not include extension.
        """
        self.db_reaction = db_reaction
        self.reaction_object = self._make_rxn_block()
        # self.reaction_container = self._make_reaction_container()
        self.metadata = services.data_export.metadata.ReactionMetaData(
            db_reaction
        ).get_dict()
        self.container_name = container_name
        self.filename = filename
        self.memory_file = None
        self.file_contents = None
        self.file_hash = None
        # if self.reaction_object:
        #     self.reaction_object.meta.update(self.metadata)

    def save(self):
        """Calls the appropriate method to save the data. Can't use .RDF without a reaction_container object"""
        if self.reaction_object:
            self.filename += ".rdf"
            self._save_as_rdf()
        else:
            self.filename += ".json"
            self._save_as_json()
        self._save_blob()

    def _make_rxn_block(self) -> Optional[str]:
        """
        Makes a rxn block with RDKIT for the
        Makes a reaction container by using the SMILES obtained from the db reaction object.
        Full reaction SMILES including reactants, agents (aka reagents and solvents), and products
        """
        reaction_smiles = self._make_reaction_smiles()
        # return None if we cannot make a reaction_container from smiles - likely to be due to no smiles present.
        if reaction_smiles:
            rxn = AllChem.ReactionFromSmarts(reaction_smiles)
            return AllChem.ReactionToRxnBlock(rxn, separateAgents=True)
        return None

    def _make_reaction_smiles(self) -> Optional[str]:
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
        if reactant_smiles and product_smiles:
            return (
                ".".join(reactant_smiles)
                + ">"
                + ".".join(reagent_smiles)
                + ".".join(solvent_smiles)
                + ">"
                + ".".join(product_smiles)
            )
        return None

    def _save_blob(self):
        """Saves the blob to Azure blob service"""
        blob_client = services.file_attachments.get_blob_client(
            self.container_name, self.filename
        )
        # Upload the blob data
        upload = io.BytesIO(self.file_contents)

        try:  # todo remove before deployment because this should only happen when debugging
            blob_client.upload_blob(upload, blob_type="BlockBlob")
        except azure.core.exceptions.ResourceExistsError:
            print("duplicate resource name error. skipping second one")
            pass

        # confirm upload
        if not blob_client.exists():
            print(f"blob {self.filename} upload failed")
            abort(401)

    def _make_rdf_contents(self) -> str:
        """Makes a .RDF reaction data file contents from a RXN block and metadata"""
        current_time = datetime.datetime.now(pytz.timezone("Europe/London")).replace(
            tzinfo=None
        )
        # add rdf header
        rdf_header = f"$RDFILE 1\n$DATM\t{current_time}\n$RFMT\n"
        # Make a multiline string from the metadata with each key-value pair formatted as $DTYPE and $DATUM
        rdf_metadata = [
            f"$DTYPE {key}\n$DATUM {value}" for key, value in self.metadata.items()
        ]
        formatted_metadata = "\n".join(rdf_metadata)
        return rdf_header + self.reaction_object + formatted_metadata

    def _save_as_rdf(self):
        """Saves the file contents as a bytearray to self.memory_file to be uploaded to Azure"""
        rdf_contents = self._make_rdf_contents()
        self.memory_file = io.StringIO()
        with self.memory_file as f:
            f.write(rdf_contents)
            # add dict / extra data
            self.file_contents = bytearray(f.getvalue(), "utf-8")

    def _save_as_json(self):
        """Saves a JSON as a bytearray to self.memory_file. This is used over rdf when no reaction smiles are present"""
        self.memory_file = io.StringIO()

        with self.memory_file as f:
            json.dump(self.metadata, f)
            self.file_contents = bytearray(f.getvalue(), "utf-8")
