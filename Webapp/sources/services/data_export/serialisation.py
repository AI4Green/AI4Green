"""For serialisation export methods into JSON format or RDF format which has Mol Blocks and serialised metadata"""
import datetime
import io
import json
import os
from typing import Optional

import pytz
import sources.services.data_export.export
from rdkit.Chem import AllChem
from sources import models, services

from . import utils


class SerialisationExport:
    """Parent class for serialisation export formats"""

    def __init__(
        self, db_reaction: models.Reaction, filename: str, container_name: str
    ):
        """
        Creates a serialisation export format and sets up the variables common across all of them.

        Args:
             db_reaction - the database entry for the reaction
             filename - the filename for the blob without the file extension. Normally the reaction_id e.g., DW1-001
             container_name - the name of the export container.
        """
        self.db_reaction = db_reaction
        self.container_name = container_name
        self.filename = filename
        self.metadata = services.data_export.metadata.ReactionMetaData(
            db_reaction
        ).get_dict()
        self.memory_file = None
        self.file_contents = None
        self.file_hash = None
        self.mime_type = None
        self.content_size = None


class RXNFileExport(SerialisationExport):
    """For exporting files as .rxn files"""

    def __init__(
        self, db_reaction: models.Reaction, filename: str, container_name: str
    ):
        """
        Creates an instance of the ReactionDataFileExport class

        Args:
             db_reaction - the database entry for the reaction
             filename - the filename for the blob without the file extension. Normally the reaction_id e.g., DW1-001
             container_name - the name of the export container.
        """
        super().__init__(db_reaction, filename, container_name)
        # Use RXN file saved in database, or generate if not found
        if self.db_reaction.reaction_rxn:
            self.reaction_object = self.db_reaction.reaction_rxn
        else:
            self.reaction_object = self._make_rxn_block()

    def save(self, extension=True):
        """Calls the appropriate method to save the data. Can't use .rxn without a reaction object"""
        if self.reaction_object:
            if extension is True:
                self.filename += ".rxn"
            self._save_as_rxn()
            self.mime_type = "chemical/x-mdl-rxn"

            self.file_hash = services.file_attachments.sha256_from_file_contents(
                self.file_contents
            )
            sources.services.data_export.export.save_blob(
                self.container_name, self.filename, self.file_contents
            )

    def _make_rxn_block(self) -> Optional[str]:
        """
        Makes a rxn block with RDKIT by using the SMILES obtained from the db reaction object.
        Full reaction SMILES including reactants, agents (aka reagents and solvents), and products
        Returns:
            RXN Block as a string
        """
        reaction_smiles = self._make_reaction_smiles()
        # return None if we cannot make a reaction_container from smiles - likely to be due to no smiles present.
        try:
            if reaction_smiles:
                # try with full reaction smiles including agents
                rxn = AllChem.ReactionFromSmarts(reaction_smiles, useSmiles=True)
            else:
                raise ValueError  # Use just reactants and products if agents fail
        except ValueError:
            try:
                rxn = AllChem.ReactionFromSmarts(
                    self.db_reaction.reaction_smiles, useSmiles=True
                )
            except ValueError:
                return None
        rxn_block = AllChem.ReactionToRxnBlock(rxn, separateAgents=True)
        return rxn_block

    def _make_reaction_smiles(self) -> Optional[str]:
        """
        Makes a reaction smiles in format reactant1.reactants2>reagents.solvents>products from a db reaction object

        Returns:
            the reaction smiles string.
        """
        reactant_smiles = utils.remove_default_data(self.db_reaction.reactants)
        reagent_smiles = utils.remove_default_data(self.db_reaction.reagents)
        solvent_smiles = services.all_compounds.get_smiles_list(
            primary_key_ls=self.db_reaction.solvent,
            person=self.db_reaction.creator_person,
        )
        product_smiles = utils.remove_default_data(self.db_reaction.products)
        if reactant_smiles and product_smiles:
            # if no reagents or solvents then the dot must be removed
            return (
                (
                    ".".join(reactant_smiles)
                    + ">"
                    + ".".join(reagent_smiles)
                    + "."
                    + ".".join(solvent_smiles)
                    + ">"
                    + ".".join(product_smiles)
                )
                .replace(">.", ">")
                .replace(".>", ">")
            )
        return None

    def _save_as_rxn(self):
        """Saves the file contents as a bytearray to self.memory_file to be uploaded to Azure"""
        self.memory_file = io.StringIO()
        with self.memory_file as f:
            f.write(self.reaction_object)
            self.file_contents = bytearray(f.getvalue(), "utf-8")
            self.content_size = f.seek(0, os.SEEK_END)


class ReactionDataFileExport(SerialisationExport):
    """For exporting files as reaction data file (.RDF) Consisting of RXN and MOL blocks."""

    # needed for reading files in, currently only used in testing but in future will want to implement ourselves
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
        Creates an instance of the ReactionDataFileExport class

        Args:
             db_reaction - the database entry for the reaction
             filename - the filename for the blob without the file extension. Normally the reaction_id e.g., DW1-001
             container_name - the name of the export container.
        """
        super().__init__(db_reaction, filename, container_name)
        self.reaction_object = RXNFileExport(
            self.db_reaction, None, None
        ).reaction_object
        self.metadata = services.data_export.metadata.ReactionMetaData(
            self.db_reaction
        ).get_dict()

    def save(self, extension=True):
        """Calls the appropriate method to save the data. Can't use .RDF without a reaction_container object"""
        if self.reaction_object:
            if extension is True:
                self.filename += ".rdf"
            self._save_as_rdf()
            self.mime_type = "chemical/x-mdl-rdfile"
        else:
            if extension is True:
                self.filename += ".json"
            self._save_as_json()
            self.mime_type = "application/json"

        self.file_hash = services.file_attachments.sha256_from_file_contents(
            self.file_contents
        )
        sources.services.data_export.export.save_blob(
            self.container_name, self.filename, self.file_contents
        )

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
            self.content_size = f.seek(0, os.SEEK_END)

    def _save_as_json(self):
        """
        Saves a JSON as a bytearray to self.memory_file.
        This is used instead of rdf if there is an error with the reaction SMILES
        """
        self.memory_file = io.StringIO()
        with self.memory_file as f:
            json.dump(self.metadata, f)
            self.file_contents = bytearray(f.getvalue(), "utf-8")
            self.content_size = f.seek(0, os.SEEK_END)


class JsonExport(SerialisationExport):
    """Class for exporting files as .JSON."""

    def __init__(
        self, db_reaction: models.Reaction, filename: str, container_name: str
    ):
        """
        Creates an instance of the JsonExport class

        Args:
             db_reaction - the database entry for the reaction
             filename - the filename for the blob without the file extension. Normally the reaction_id e.g., DW1-001
             container_name - the name of the export container.
        """
        super().__init__(db_reaction, filename, container_name)

    def save(self, extension=True):
        self._save_json()
        if extension is True:
            self.filename += ".json"
        self.mime_type = "application/json"
        self.file_hash = services.file_attachments.sha256_from_file_contents(
            self.file_contents
        )
        sources.services.data_export.export.save_blob(
            self.container_name, self.filename, self.file_contents
        )

    def _save_json(self):
        """Saves a JSON as a bytearray to self.memory_file."""
        self.memory_file = io.StringIO()

        with self.memory_file as f:
            json.dump(self.metadata, f)
            self.file_contents = bytearray(f.getvalue(), "utf-8")
            self.content_size = f.seek(0, os.SEEK_END)
