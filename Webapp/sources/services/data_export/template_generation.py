import io
import json
import re
from datetime import datetime
from docx import Document
from flask import render_template
from io import StringIO
from sources import models, services

from . import utils


class PdfExport:
    """Class for exporting files as .PDF."""

    def __init__(
        self, db_reaction: models.Reaction, filename: str, container_name: str
    ):
        """
        Creates an instance of the PdfExport class

        Args:
             db_reaction - the database entry for the reaction
             filename - the filename for the blob without the file extension. Normally the reaction_id e.g., DW1-001
             container_name - the name of the export container.
        """
        self.db_reaction = db_reaction
        self.metadata = services.data_export.metadata.ReactionMetaData(
            db_reaction
        ).get_dict()
        self.container_name = container_name
        self.filename = filename
        self.memory_file = None
        self.file_contents = None
        self.file_hash = None

    def make(self):
        """Generates a PDF from the summary html template"""
        summary_data = json.loads(self.db_reaction.summary_table_data)
        reaction_table_data = json.loads(self.db_reaction.reaction_table_data)
        reaction_data = {**summary_data, **reaction_table_data}
        # remove empty values
        reaction_data = process_reaction_data(reaction_data)
        # physical form num to str needed

        reactant_data = services.summary.get_reactant_data(reaction_data)
        print(reactant_data)
        document = Document()
        document.add_heading("Sample Press Release", 0)
        f = StringIO()
        document.save(f)
        f.seek(0)

    def save(self):
        self._save_json()
        self.filename += ".pdf"
        services.data_export.export.save_blob(
            self.container_name, self.filename, self.file_contents
        )

    def _save_json(self):
        """Saves a JSON as a bytearray to self.memory_file. This is used over rdf when no reaction smiles are present"""
        self.memory_file = io.StringIO()

        with self.memory_file as f:
            json.dump(self.metadata, f)
            self.file_contents = bytearray(f.getvalue(), "utf-8")


def process_reaction_data(reaction_data):
    reaction_data = {
        key.replace("_names", "s"): value for key, value in reaction_data.items()
    }
    reaction_data = {
        update_rounded_raw_keys(key): value for key, value in reaction_data.items()
    }
    reaction_data = process_physical_forms_keys(reaction_data)
    # conditional changes
    for key, value in reaction_data.items():
        if isinstance(value, list):
            reaction_data[key] = remove_trailing_empty_strings(value)
    return reaction_data


def process_physical_forms_keys(data):
    keys_to_remove = []
    for key, value in data.items():
        if key.endswith("_physical_forms_text"):
            base_key = key[:-5]  # Remove '_text' from the key
            # if base_key in data:
            #     del data[base_key]  # Remove corresponding '_physical_forms' key
            # data[base_key] = value.rstrip('_text')  # Update value without '_text'
            data[base_key] = value
            keys_to_remove.append(key)
    for key in keys_to_remove:
        if key in data:
            del data[key]
    return data


def update_rounded_raw_keys(string: str) -> str:
    # Pattern to match *_masses, *_amounts, *_volumes, etc.
    rounded_pattern = re.compile(r"^(.*)_(masses|amounts|volumes|)$")
    raw_pattern = re.compile(r"^(.*)_(masses|amounts|volumes|)_raw$")

    # Check if the string matches *_raw pattern
    raw_match = raw_pattern.match(string)
    if raw_match:
        return f"{raw_match.group(1)}_{raw_match.group(2)}"

    # Check if the string matches the base pattern
    rounded_match = rounded_pattern.match(string)
    if rounded_match:
        return f"rounded_{string}"

    return string


def remove_trailing_empty_strings(lst: list) -> list:
    while lst and lst[-1] == "":
        lst.pop()
    return lst
