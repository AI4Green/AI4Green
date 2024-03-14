import ast
import json

import chython.files
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

    def __init__(self, db_reaction: models.Reaction, filename: str):
        """
        Filename should not include extension.
        """
        self.db_reaction = db_reaction
        self.reaction_container = self.make_reaction_container()
        self.metadata = services.data_export.metadata.ReactionMetaData(
            db_reaction
        ).get_dict()
        self.filename = filename
        if self.reaction_container:
            self.reaction_container.meta.update(self.metadata)

    def make_reaction_container(self):
        """
        Makes a reaction container by using the SMILES obtained from the db reaction object.
        Full reaction SMILES including reactants, agents (aka reagents and solvents), and products
        """
        reaction_smiles = self.make_reaction_smiles()
        # return None if we cannot make a reaction_container from smiles - likely to be due to no smiles present.
        try:
            return chython.files.smiles(reaction_smiles)
        except ValueError:
            return None

    def make_reaction_smiles(self) -> str:
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

    def save(self):
        """Calls the appropriate method to save the data. Can't use .RDF without a reaction_container object"""
        if self.reaction_container:
            self.filename += ".rdf"
            self.save_as_rdf()
        else:
            self.filename += ".json"
            self.save_as_json()

    def save_as_rdf(self):
        """Saves as an RDF"""
        with chython.files.RDFWrite(self.filename) as f:
            f.write(self.reaction_container)

    def literal_eval_metadata(self, rdf_contents: chython.ReactionContainer):
        """Read the rdf values literally to reintroduce their types"""

        for key in self.metadata.keys():
            # if the string is equal to none we reload
            if (
                rdf_contents.meta[key] == "None"
                or key in self.non_string_types_inside_meta
                and isinstance(rdf_contents.meta[key], str)
            ):
                rdf_contents.meta.update(
                    {key: ast.literal_eval(rdf_contents.meta[key])}
                )

    def save_as_json(self):
        with open(self.filename, "w") as f:
            json.dump(self.metadata, f)

        with open(self.filename, "r") as f:
            meta_reload = json.load(f)

        assert meta_reload == self.metadata, "change in data during file read/write"
