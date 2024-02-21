import json
import re
from typing import Dict, List, Tuple

import chython.files

# from chython import files as ctf  # - chemical table files
from chython import smiles
from chython.containers import MoleculeContainer, ReactionContainer
from flask import abort, jsonify, request
from rdkit.Chem import AllChem
from sources import auxiliary, db, models, services


def make_reaction_smiles(reaction: models.Reaction) -> str:
    """
    Uses an instance of a reaction object form the database to make
     a reaction smiles string in format reactant1.reactants2>reagents.solvents>products

    Args:
        reaction - the reaction from which we are making a reaction smiles string

    Returns:
        the reaction smiles string.
    """
    reactant_smiles = reaction.reactants
    reagent_smiles = services.all_compounds.get_smiles_list(reaction.reagents)
    solvent_smiles = services.all_compounds.get_smiles_list(reaction.solvent)
    product_smiles = reaction.products
    reaction_smiles = (
        ".".join(reactant_smiles)
        + ">"
        + ".".join(reagent_smiles)
        + ".".join(solvent_smiles)
        + ">"
        + ".".join(product_smiles)
    )
    return reaction_smiles


#
# def add_metadata(
#     db_reaction: models.Reaction,
#     reaction_container: chython.containers.ReactionContainer,
# ):
#     """
#     Adds metadata to the reaction_container for inclusion in the .reaction-data-file (.rdf) file export for a reaction.
#     We add this by updating the 'meta' dictionary property of the chython ReactionContainer object
#
#     Args:
#         db_reaction: The database object for the reaction
#         reaction_container: The chython ReactionContainer object which can be written to a .rdf file type.
#     """
#     rdf_metadata = [
#         "temperature",
#         "solvents",
#         "reagents",
#         "creator_username",
#         "creator_workbook",
#         "creator_workgroup",
#         "time_of_creation",
#         "reaction_completed",
#         "purification_method",
#         "batch_or_flow",
#         "percent_yield",
#         "reactants",
#         "experimental_writeup",
#         "standard_protocols_used",
#         "nmr_data",
#         "file_attachment_names",
#     ]
#
#     rxn_data = json.loads(db_reaction.reaction_table_data)
#     summary_data = json.loads(db_reaction.summary_table_data)
#     solvent_names = (
#         [x for x in rxn_data["solvent_names"]] if rxn_data.get("solvents") else None
#     )
#     solvent_smiles = services.all_compounds.get_smiles_list(
#         rxn_data.get("solvent_primary_keys")
#     )
#     # solvent_smiles = [x for x in solvent_smiles] if rxn_data.get('solvents') else None
#     solvent_volumes = (
#         [x for x in rxn_data["solvent_volumes"]] if rxn_data.get("solvents") else None
#     )
#
#     # solvent_dict = [{'name': name, 'SMILES': smiles, 'volume': volume, 'unit': solvent_unit} for name, smiles, volume in solvent_names, solvent_smiles, solvent_volume]
#     # reagent_dict = {'name': name, 'SMILES': smiles, 'amount': amount, 'unit': unit}
#     # #
#     # solvent_dict = {'solvents': {'name':}}
#
#     update_dict = {
#         "temperature": db_reaction.summary_table_data.get("temperature", None),
#         "solvents": db_reaction.reaction_table_data.get("solvents", None),
#     }
#
#     reaction_container.meta.update()


# TODO
def make_rxn_file(reaction: models.Reaction):
    # reaction_smiles = "F[B-](F)(F)F.[Zn+]C1=CC=CC=C1.[Zn](C1=CC=CC=C1)C1=CC=CC=C1>>O |f:0.1,c:7,9,14,16,21,23,t:5,12,19,lp:0:3,2:3,3:3,4:3,25:2|"
    # reaction = services.reaction.get_from_reaction_id_and_workbook_id('DW1-001', 1)
    dbreaction = (
        db.session.query(models.Reaction).filter(models.Reaction.id == 1).first()
    )

    reaction_smiles = make_reaction_smiles(dbreaction)
    print(reaction_smiles)
    reaction = chython.files.smiles(reaction_smiles)
    # add_metadata(dbreaction, reaction)
    reaction.meta.update({"temp": "hot"})

    # rxn_contents = AllChem.ReactionToRxnBlock(rxn, separateAgents=True)
    dummy_file = "dummy_rxn.rdf"

    with chython.files.RDFWrite(dummy_file) as f:  # context manager supported
        # for r in data:
        f.write(reaction)

    with chython.files.RDFRead(dummy_file) as f:
        first = next(f)
        print(first)

        # # create a RDFParser object, this is a generator that yields Reaction objects
        # rdfreader = RDFParser(
        #     rdf_file,
        #     # except_on_invalid_molecule=False,
        #     # will return None instead of raising an exception if a molecule is invalid
        #     # except_on_invalid_reaction=False,
        #     # will return None instead of raising an exception if a reaction is invalid
        # )
        #
        # for rxn5 in rdfreader:
        #     print(rxn5)
        #     # rxn is a Reaction object, it is several attributes, including:
        #     print(rxn5.smiles)  # reaction SMILES string
        #     print(rxn5.properties)  # a dictionary of properties extracted from the RXN record
        #
        #     reactants1 = rxn5.reactants  # a list of Molecule objects
        #     products1 = rxn5.products
        #     solvents1 = rxn5.solvents
        #     catalysts1 = rxn5.catalysts
        #
        #     # Molecule objects have several attributes, including:
        #     print(reactants1[0].smiles)
        #     print(reactants1[0].properties)  # a dictionary of properties extracted from the MOL record (often empty)
        #     print(reactants1[0].rd_mol)  # an RDKit molecule object

    # # to test
    # with open("dummy_rxn.rxn", "r") as file:
    #     file_contents = file.read()
    #     rxn2 = AllChem.ReactionFromRxnBlock(file_contents)
    #
    # AllChem

    # test_file = "dummy_rxn.rdf"
    # rxn2 = AllChem.ReactionFromRxnFile(test_file)
    #
    # rxn1smiles = AllChem.ReactionToSmiles(rxn)
    # rxn2smiles = AllChem.ReactionToSmiles(rxn2)
    #
    # print(rxn1smiles, "\n", rxn2smiles)
    #
    # if rxn == rxn2:
    #     print("yay :)")
    # else:
    #     print("nay :(")

    # return {}
