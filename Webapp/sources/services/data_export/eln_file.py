import io
import json
from datetime import datetime
from itertools import zip_longest
from typing import Dict, Generator, List, Optional, Tuple

import pytz
import sources.services.data_export.export
from flask import request
from rdkit.Chem import AllChem
from sources import models, services

from . import utils

#     """
#     We are making a .eln file for a specific workbook. This is a zipped directory containing a ro-crate-metadata.json
#     to describe the contents.
#     """
#
#     # reaction_list = services.reaction.list_active_in_workbook(
#     #     workbook, workgroup, sort_crit="time"
#     # )
#     workgroup_name = request.form["workgroup"]
#     workbook_name = request.form["workbook"]
#     eln_file = ELNFile(workgroup_name, workbook_name)
#     # eln_file = DataExport(workgroup_name, workbook_name).to_eln_file()
#     return eln_file

# made according to this specification https://github.com/TheELNConsortium/TheELNFileFormat
# first we describe the ro-crate metadata_json
# ro_crate_metadata_json_contents = describe_ro_crate_metadata_json()
# now for each reaction we want to make a research object crate
# for idx, reaction in enumerate(reaction_list):
# make a folder per experiment.


# def get_existing_files(self):
#     file_attachments =


class ELNFileExport:
    # absolute basic just include an rdf for a reaction

    def __init__(self, data_export_request: models.DataExportRequest):
        self.data_export_request = data_export_request
        self.reaction_list = None
        self.time = datetime.now(pytz.timezone("Europe/London")).replace(tzinfo=None)
        self.parts = []

    def make_eln_file(self):
        # name the root?
        root = (
            self.data_export_request.workgroup.name
            + "-"
            + self.data_export_request.workbooks[0].name
            + "-"
            + self.data_export_request.time_of_request.strftime("%Y-%m-%d-%H%M%S")
            + "-export"
        )
        print(root)
        # for each reaction extract and name the parts needed to make the ELN file

        # define authors
        self._define_authors()
        for reaction in self.data_export_request.reactions:
            export_rxn = ELNExportReaction(reaction, self.data_export_request)
            export_rxn.get_files()
            export_rxn.get_parts()
            self.reaction_list.append(export_rxn)
            # get the parts
            self.parts.append()

        self._make_ro_crate_metadata_json()

    def _make_ro_crate_metadata_json(self):
        version, git_hash = services.utils.get_app_version()
        # this part is common to all .eln files from AI4Green
        metadata = {
            "@context": "https://w3id.org/ro/crate/1.1/context",
            "@graph": [
                {
                    "@id": "ro-crate-metadata.json",
                    "@type": "CreativeWork",
                    "about": {"@id": "./"},
                    "conformsTo": {"@id": "https://w3id.org/ro/crate/1.1"},
                    "dateCreated": datetime.now(pytz.timezone("Europe/London")).replace(
                        tzinfo=None
                    ),
                    "sdPublisher": {
                        "@type": "Organization",
                        "areaServed": "Pisces–Cetus Supercluster Complex",
                        "name": "AI4Green",
                        # "logo": "https://www.elabftw.net/img/elabftw-logo-only.svg",
                        "slogan": "AI for Green Chemistry. Electronic Lab Notebook with Machine Learning Support.",
                        "url": "https://www.ai4green.app",
                        "parentOrganization": {
                            "@type": "Organization",
                            "name": "University of Nottingham",
                            "logo": "https://www.nottingham.ac.uk/Brand/LegacyAssets/images-multimedia/2022/Logos/BrandEvolution-NottinghamBlue-Cropped-450x173.png",
                            # "slogan": "Open Source software for research labs.",
                            "url": "https://www.nottingham.ac.uk/",
                        },
                    },
                    "version": "1.0",
                },
                {
                    "@id": "#ro-crate_created",
                    "@type": "CreateAction",
                    "object": {"@id": "./"},
                    "name": "RO-Crate created",
                    "endTime": datetime.now(pytz.timezone("Europe/London")).replace(
                        tzinfo=None
                    ),
                    "instrument": {
                        "@id": "https://www.ai4green.app",
                        "@type": "SoftwareApplication",
                        "name": "AI4Green",
                        "version": version,
                        "git_commit_hash": git_hash,
                        "identifier": "https://www.ai4green.app",
                    },
                    "actionStatus": {"@id": "http://schema.org/CompletedActionStatus"},
                },
                # now from here this is export specific.
                self._describe_root_directory(),
                # for ELNExportFile in ELNExportFiles todo
                {
                    "@id": "./2025-01-02 - Synthesis-of-Aspirin - 64e381cf/export-elabftw.json",
                    "@type": "File",
                    "description": "JSON export",
                    "name": "export-elabftw.json",
                    "encodingFormat": "application/json",
                    "contentSize": "5137",
                    "sha256": "75dc6abdeba145b4c7586ed2df4145427cadec992c34dae52dde927fba972c8a",
                },
                # for ELNExportDataset in ELNExportDatasets
                {
                    "@id": "./2025-01-02 - Synthesis-of-Aspirin - 64e381cf",
                    "@type": "Dataset",
                    "author": {
                        "@id": "person://9d4b405329eb4659017b896f33c8501327e7acbc11ef18b325f7b18770f00dde?hash_algo=sha256"
                    },
                    "dateCreated": "2024-04-30T02:23:27+02:00",
                    "dateModified": "2024-04-30T11:10:04+02:00",
                    "identifier": "20240430-64e381cfd2039556c243efafe146e9ab867015fc",
                    "comment": [],
                    "keywords": "chemistry,has-mathjax",
                    "name": "Synthesis of Aspirin",
                    "text": "Hello Pisces–Cetus Supercluster Complex",
                    "url": "https://demo.elabftw.net/experiments.php?mode=view&id=262",
                    # Contents of dataset
                    "hasPart": [
                        {
                            "@id": "./2025-01-02 - Synthesis-of-Aspirin - 64e381cf/export-elabftw.json"
                        }
                    ],
                    "mentions": [],
                    "category": "Production",
                    "status": "Success",
                },
                # Contents of ro-crate
                {
                    "@id": "./",
                    "@type": "Dataset",
                    "hasPart": [
                        {"@id": "./2025-01-02 - Synthesis-of-Aspirin - 64e381cf"}
                    ],
                },
                # define comments
                {
                    "@id": "comment://2023-09-28T03%3A44%3A52%2B02%3A00",
                    "@type": "Comment",
                    "dateCreated": "2023-09-28T03:44:52+02:00",
                    "text": "Well, it&#39;s a relief to see that the scientific community is hard at work discovering what we already suspected: that the sky is indeed blue. I can&#39;t wait for their next groundbreaking revelation that grass is green.",
                    "author": {
                        "@id": "person://a0568e4aa35eec2cf3962db57eec039a205d3b7d7c7132388af976dfaf8f2694?hash_algo=sha256"
                    },
                },
                # define authors
                {
                    "@id": "person://9d4b405329eb4659017b896f33c8501327e7acbc11ef18b325f7b18770f00dde?hash_algo=sha256",
                    "@type": "Person",
                    "familyName": "Wiegand",
                    "givenName": "Katelyn",
                },
            ],
        }
        print(metadata)

        # make the ro crate metadata json first todo - update with list of reaction folders

        # for each reaction make a folder

    def _describe_root_directory(self):
        return {
            "@id": "./",
            "@type": ["Dataset"],
            "hasPart": [part.id for part in self.parts],
        }

    def get_author(self, reaction: models.Reaction) -> Dict:
        person = reaction.creator_person
        if person.user:
            author = {
                "@id": f"./author/{person.id}",
                "@type": "Person",
                "name": person.user.fullname,
                "email": person.user.email,
                "workgroup_role": person.user,
            }
        else:
            author = {
                "@id": f"./author/{person.id}",
                "@type": "Person",
                "name": "deleted",
                "email": "deleted",
                "workgroup_role": person.user,
            }
        return author

    def describe_ro_crate_metadata(self) -> Dict:
        """
        Returns the contents to make the metadata json file.
        """
        return {}

    def get_has_parts_id_list(reaction: models.Reaction) -> List[Dict]:
        """
        Returns a list of dictionaries, 1 per hasPart / file attachments and other files generated for export.
        """
        has_parts_list = []
        for file in reaction.file_attachments:
            has_parts_list.append({"@id": f"file-id/{reaction.reaction_id}/{file.id}"})
        # include any MOL, SDF or RDF files. - id would follow same pattern but swap file.id for file extension.
        return has_parts_list

    def get_file_details(file_id: str, reaction: models.Reaction) -> Tuple[Dict, str]:
        """
        Returns a string describing the file from its extension type
        """
        file_details = {}
        name = ""
        return file_details, name

    def get_has_parts_definitions(
        self, reaction: models.Reaction, id_list: List[Dict]
    ) -> Generator:
        """
        Returns a Tuple of dictionaries containing the definitions of the hasParts which were identified earlier
        """
        definition_list = []
        for file_id, file in zip_longest(id_list, reaction.file_attachments):
            if file:
                file_details = file.file_details
                name = file.name
            else:
                name, file_details = self.get_file_details(file_id["@id"], reaction)

            definition_list.append(
                {
                    "@id": file_id["@id"],
                    "@type": "File",
                    "description": "",  # Optional?
                    "name": name,
                    "encodingFormat": file_details["mimetype"],
                    "contentSize": file_details["size"],
                    "sha256": file_details[
                        "sha256"
                    ],  # TODO - need to make process to generate this somewhere.
                }
            )
        return (definition for definition in definition_list)

        #  also remember one of: [ .mol, .SDF, .RDF], csv with headings required to fill in the table, experimental.txt, pdf summary of rxn if not already


class ELNExportReaction:
    def __init__(self, reaction: models.Reaction, eln_export: ELNFileExport):
        self.reaction = reaction
        self.data_export_request = eln_export.data_export_request
        self.time = eln_export.time
        self.files = []
        self.defined_comments = []
        self.defined_files = []
        self.parts = []

    def get_files(self):
        # for now we JUST make rdf
        rdf = services.data_export.serialisation.ReactionDataFileExport(
            self.reaction, self.reaction.reaction_id, self.data_export_request.uuid
        )
        rdf.save()
        rdf_dict = {
            "filetype": "RDF",
            "made_for_export": True,
            "filename": rdf.filename,
            "container": rdf.container_name,
            "size": "BIG",
            "sha256": rdf.file_hash,
            "mimetype": rdf.mime_type,
        }
        self.files.append(rdf_dict)

        # get pre existing files
        pre_existing_files = self.reaction.file_attachments.to_dict()
        for file in pre_existing_files:
            file["made_for_export"] = False

        self.files.append(pre_existing_files)

    def get_parts(self):
        """Make all parts for a reaction including the dataset and the files"""
        self._define_comments()
        self._define_files()
        self._define_reaction_dataset()

    def _define_reaction_dataset(self):
        # each reaction has 1 part for the dataset (reaction) and n number of files for the reaction
        self.parts.append(
            {
                "@id": self.reaction.reaction_id,
                "@type": "Dataset",
                "name": self.reaction.name,
                # "identifier": self.reaction.uuid, reaction has no uuid yet.
                "author": {"@id": f"./author/{self.reaction.creator_person.id}"},
                "dateCreated": self.reaction.time_of_creation,
                "dateModified": self.reaction.time_of_update,
                "comment": self._get_comment_id_list(),
                "hasPart": self._get_file_part_ids(),  # files
                "url": f"{request.base_url}/{self.reaction.workgroup}/{self.reaction.workbook}/{self.reaction.reaction_id}/no",
                # may remove no in future
            }
        )

    def _get_file_part_ids(self) -> List[Dict]:
        """
        Returns a list of dictionaries, 1 per hasPart / file attachments and other files generated for export.
        """
        has_parts_list = []
        for file in self.files:
            has_parts_list.append(
                {"@id": f"file-id/{self.reaction.reaction_id}/{file['id']}"}
            )
        # include any MOL, SDF or RDF files. - id would follow same pattern but swap file.id for file extension.
        return has_parts_list

    def _define_comments(self):
        for comment in self.reaction.addenda:
            self.defined_comments.append(
                {
                    "@id": f"comment-id/{self.reaction.reaction_id}/{comment.id}",
                    "@type": "Comment",
                    "dateCreated": comment.time_of_creation,
                    "text": comment.text,
                    "author": {"@id": f"./author/{self.reaction.creator_person.id}"},
                }
            )

    def _define_files(self):
        for file in self.files:
            self.defined_files.append(
                {
                    "@id": "./IntermetalsAtInterfaces/001_GetSprayMachine/metadata.json",
                    "@type": "File",
                    "name": "export_metadata.json",
                    "description": "JSON export",
                    "encodingFormat": "application/json",
                    "contentSize": "338",
                    "sha256": "0aa2d3acd0d6bfab8d05e33fc45824d05e02eb1863f94d1cb8cd8471180e7947",
                    "dateModified": "2024-02-18T15:49:00.325602",
                }
            )

    # def _get_comment_id_list(self) -> List[Dict]:
    #     """
    #     Returns a list of dictionaries, 1 per comment containing the id of the comment
    #     """
    #     comment_id_list = []
    #     for comment in self.reaction.addenda:
    #         # make an id from the work
    #         comment_id_list.append(
    #             {"@id": f"comment-id/{self.reaction.reaction_id}/{comment.id}"}
    #         )
    #     return comment_id_list

    # def make_files(self, reaction):
    #     """
    #     Saves files to tmp storage
    #     """
    # make .rdf

    # file json example

    # {
    #   "@id": "./2025-03-03 - Testing-the-eLabFTW-lab-notebook - b6db898e/example.png",
    #   "@type": "File",
    #   "description": "",
    #   "name": "example.png",
    #   "alternateName": "6b/6befead611faa1c99829ee811f5973e308f28cc8faea98594efc9ecf637bb62db3da3172e5c53e23f03bec7750b7ead5a32542b099544f8ee2902e4d66a99583.png",
    #   "contentSize": 40959,
    #   "sha256": "f2780bb5a8883f8c89a6da1c5bf4604f7cb44bdcb093e16484eb81c1bcb1f580"
    # },
    #
    #
    #
    #
    #
    # {
    #     "@id": "./2025-03-03 - Testing-the-eLabFTW-lab-notebook - b6db898e",
    #     "@type": "Dataset",
    #     "author": {
    #         "@id": "person://a0568e4aa35eec2cf3962db57eec039a205d3b7d7c7132388af976dfaf8f2694?hash_algo=sha256"
    #     },
    #     "dateCreated": "2023-09-28T03:44:54+02:00",
    #     "dateModified": "2023-10-02T14:51:49+02:00",
    #     "identifier": "20230928-b6db898e231dbc2234227d19656d7e19a77654b7",
    #     "comment": [
    #         {
    #             "@id": "comment://2023-09-28T03%3A44%3A54%2B02%3A00"
    #         }
    #     ],
    #     "keywords": [
    #         "demo",
    #         "test"
    #     ],
    #     "name": "Testing the eLabFTW lab notebook",
    #     "text": "<h1>Goal</h1>\n<p>Test the software.</p>\n<h1>Procedure</h1>\n<p>Click everywhere and explore everything.</p>\n<h1>Results</h1>\n<p>It's really nice, I think I'll adopt it for our lab.</p>",
    #     "url": "https://elab.local:3148/experiments.php?mode=view&id=264",
    #     "hasPart": [
    #         {
    #             "@id": "./2025-03-03 - Testing-the-eLabFTW-lab-notebook - b6db898e/export-elabftw.json"
    #         },
    #         {
    #             "@id": "./2025-03-03 - Testing-the-eLabFTW-lab-notebook - b6db898e/example.png"
    #         },
    #         {
    #             "@id": "./2025-03-03 - Testing-the-eLabFTW-lab-notebook - b6db898e/example.pdf"
    #         }
    #     ],
    #     "mentions": [
    #         {
    #             "@id": "https://elab.local:3148/database.php?mode=view&id=51"
    #         },
    #         {
    #             "@id": "https://elab.local:3148/database.php?mode=view&id=50"
    #         },
    #         {
    #             "@id": "https://elab.local:3148/experiments.php?mode=view&id=256"
    #         },
    #         {
    #             "@id": "https://elab.local:3148/experiments.php?mode=view&id=257"
    #         },
    #         {
    #             "@id": "https://elab.local:3148/experiments.php?mode=view&id=258"
    #         }
    #     ]
    # },

    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    # {
    #     "@id": "./some-unique-id/23",
    #     "@type": "Dataset",
    #     "author": {
    #         "@id": "./some-author-id/44"
    #     },
    #     "dateCreated": "2023-09-23T01:02:26+02:00",
    #     "dateModified": "2023-09-27T23:02:44+02:00",
    #     "comment": [
    #         {
    #             "@id": "./some-comment-id/91"
    #         }
    #     ],
    # },
    # {
    #     "@id": "./some-comment-id/91",
    #     "@type": "Comment",
    #     "dateCreated": "2023-09-23T01:02:26+02:00",
    #     "text": "This is the content of the comment.",
    #     "author": {
    #         "@id": "./some-author-id/44"
    #     }
    # },
    # {
    #     "@id": "./some-author-id/44"
    #            "@type": "Person",
    # "familyName": "Tapie",
    # "givenName": "Bernard"
    # }
    # }
    #
    #
