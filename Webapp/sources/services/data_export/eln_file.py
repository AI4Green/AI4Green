import io
import json
import uuid
from datetime import datetime
from itertools import zip_longest
from typing import Dict, Generator, List, Optional, Tuple
from urllib.parse import quote

import pytz
import sources.services.data_export.export
from flask import current_app, request, url_for
from rdkit.Chem import AllChem
from sources import models, services


class ELNFileExport:
    # absolute basic just include an rdf for a reaction

    def __init__(self, data_export_request: models.DataExportRequest):
        self.data_export_request = data_export_request
        self.root = (
            # self.data_export_request.WorkGroup.name
            # + "-"
            self.data_export_request.workbooks[0].name
            + "-"
            + self.data_export_request.time_of_request.strftime("%Y-%m-%d")
            + "-export"
        )
        self.reaction_list = []
        self.time = (
            datetime.now(pytz.timezone("Europe/London"))
            .replace(tzinfo=None)
            .strftime("%Y-%m-%d %H:%M:%S")
        )
        self.root_parts = []
        self.defined_authors = []

    def make_eln_file(self):
        # name the root?

        # for each reaction extract and name the parts needed to make the ELN file

        # define authors
        self._define_authors()
        for reaction in self.data_export_request.reactions:
            export_rxn = ELNExportReaction(reaction, self)
            export_rxn.get_files()
            export_rxn.get_parts()
            self.reaction_list.append(export_rxn)
            # get the root parts
            self.root_parts = export_rxn.defined_datasets

        self._make_ro_crate_metadata_json()

    def _make_ro_crate_metadata_json(self):
        # TODO Rework this to make a list of dicts to exactly match the output

        ro_crate_graph = []
        ro_crate_graph += get_constant_ro_crate_start()
        ro_crate_graph.append(self._describe_root_directory())
        ro_crate_graph += [
            comment for rxn in self.reaction_list for comment in rxn.defined_comments
        ]
        ro_crate_graph += [
            file for rxn in self.reaction_list for file in rxn.defined_files
        ]
        ro_crate_graph += [
            dataset for rxn in self.reaction_list for dataset in rxn.defined_datasets
        ]
        ro_crate_graph += [author for author in self.defined_authors]

        # ro_crate_graph.append(get_constant_ro_crate_start())
        # ro_crate_graph.append(self._describe_root_directory())
        # ro_crate_graph.append([comment for rxn in self.reaction_list for comment in rxn.defined_comments])
        # ro_crate_graph.append([file for rxn in self.reaction_list for file in rxn.defined_files])
        # ro_crate_graph.append([dataset for rxn in self.reaction_list for dataset in rxn.defined_datasets])
        # ro_crate_graph.append([author for author in self.defined_authors])

        # defined_ro_crate_json = get_constant_ro_crate_start()
        # defined_root_dir = self._describe_root_directory()
        # defined_comments = []
        # defined_files = []
        # defined_datasets = []
        #
        # for rxn in self.reaction_list:
        #     for comment in rxn.defined_comments:
        #         defined_comments.append(comment)
        #     for file in rxn.defined_files:
        #         defined_files.append(file)
        #     for dataset in rxn.defined_datasets:
        #         defined_datasets.append(dataset)
        # # ', '.join(json.dumps(d) for d in list_of_dicts)
        # defined_authors = {}
        # [defined_authors.update(author) for author in self.defined_authors]
        # defined_authors = {}.update(author) for author in self.defined_authors)
        # this part is common to all .eln files from AI4Green
        ro_crate_metadata_contents = {
            "@context": "https://w3id.org/ro/crate/1.1/context",
            "@graph": [
                ro_crate_graph
                # defined_ro_crate_json,
                # # now from here this is export specific.
                # # contents of ro_crate
                # defined_root_dir,
                # # for ELNExportDataset in ELNExportDatasets
                # defined_datasets,
                # # for ELNExportFile in ELNExportFiles
                # defined_files,
                # # define comments
                # defined_comments,
                # # define authors
                # defined_authors,
            ],
        }
        with open("ro-crate-metadata.json", "w") as f:
            f.write(json.dumps(ro_crate_metadata_contents))

        # make the ro crate metadata json first todo - update with list of reaction folders

        # for each reaction make a folder

    def _describe_root_directory(self):
        return {
            "@id": "./",
            "@type": ["Dataset"],
            "hasPart": [{"@id": part["@id"]} for part in self.root_parts],
        }

    def _define_authors(self):
        # get the json definition of every author of a reaction in the export list
        [
            self.defined_authors.append(self.get_author(reaction))
            for reaction in self.data_export_request.reactions
        ]
        # remove duplicates
        self.defined_authors = [
            i
            for n, i in enumerate(self.defined_authors)
            if i not in self.defined_authors[:n]
        ]

    def get_author(self, reaction: models.Reaction) -> Dict:
        person = reaction.creator_person
        if person.user:
            author = {
                "@id": f"./author/{person.id}",
                "@type": "Person",
                "fullName": person.user.fullname,
                "emailAddress": person.user.email,
            }
        else:
            author = {
                "@id": f"./author/{person.id}",
                "@type": "Person",
                "fullName": "deleted",
                "emaiLAddress": "deleted",
            }
        return author


class ELNExportReaction:
    def __init__(self, reaction: models.Reaction, eln_export: ELNFileExport):
        self.reaction = reaction
        self.root = eln_export.root
        self.data_export_request = eln_export.data_export_request
        self.time = eln_export.time
        self.files = []
        self.defined_comments = []
        self.defined_files = []
        self.defined_datasets = []

    def get_files(self):
        # for now we JUST make rdf
        rdf = services.data_export.serialisation.ReactionDataFileExport(
            self.reaction, self.reaction.reaction_id, self.data_export_request.uuid
        )
        rdf.save()
        rdf_dict = {
            "filetype": "RDF",
            "made_for_export": True,
            "display_name": rdf.filename,
            "container": rdf.container_name,
            "content_size": rdf.content_size,
            "sha256": rdf.file_hash,
            "mimetype": rdf.mime_type,
            "time_of_creation": self.time,
            "uuid": uuid.uuid4(),
        }
        self.files.append(rdf_dict)

        # get pre existing files
        pre_existing_files = [
            x.to_export_dict() for x in self.reaction.file_attachments
        ]
        for file in pre_existing_files:
            file["made_for_export"] = False
        self.files += pre_existing_files

    def get_parts(self):
        """Define all parts for a reaction including the comments, files, and dataset"""
        self._define_comments()
        self._define_files()
        self._define_reaction_dataset()

    def _define_reaction_dataset(self):
        # each reaction has 1 part for the dataset (reaction) and n number of files for the reaction
        self.defined_datasets.append(
            {
                "@id": f"{self.root}/{self.reaction.reaction_id}",
                "@type": "Dataset",
                "name": self.reaction.name,
                # "identifier": self.reaction.uuid, reaction has no uuid yet.
                "author": {"@id": f"./author/{self.reaction.creator_person.id}"},
                "dateCreated": self.reaction.time_of_creation.strftime(
                    "%Y-%m-%d %H:%M:%S"
                ),
                "dateModified": self.reaction.time_of_update.strftime(
                    "%Y-%m-%d %H:%M:%S"
                ),
                "comment": [comment["@id"] for comment in self.defined_comments],
                "hasPart": [file["@id"] for file in self.defined_files],  # files
                "url": url_for("main.index", _external=True)
                + quote(
                    f"/{self.reaction.workbook.WorkGroup.name}/{self.reaction.workbook.name}/{self.reaction.reaction_id}/no"
                ),
                # may remove no in future
            }
        )

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
                    "@id": f"{self.root}/{self.reaction.reaction_id}/{file['uuid']}",
                    "@type": "File",
                    "name": f"{file['display_name']}",
                    "description": self._get_file_description(file),
                    "encodingFormat": file["mimetype"],
                    "contentSize": str(file["content_size"]),
                    "sha256": file["sha256"],
                    "dateModified": file["time_of_creation"],
                }
            )

    def _get_file_description(self, file: Dict) -> str:
        if file["made_for_export"] is True:
            export_file_type_description_dict = {
                "RDF": "A Chemistry specific format with serialized metadata and molecules represented in a "
                "Chemical Table format and assigned as reactants, products, or agents",
            }
            description = export_file_type_description_dict[file["filetype"]]
        else:
            if file["autogenerated"] is True:
                description = "Autogenerated PDF of experiment"
            else:
                description = "Reaction data uploaded by user"
        return description


def get_constant_ro_crate_start():
    version, git_hash = services.utils.get_app_version()
    return [
        {
            "@id": "ro-crate-metadata.json",
            "@type": "CreativeWork",
            "about": {"@id": "./"},
            "conformsTo": {"@id": "https://w3id.org/ro/crate/1.1"},
            "schemaVersion": "1.0",
            "dateCreated": datetime.now(pytz.timezone("Europe/London"))
            .replace(tzinfo=None)
            .strftime("%Y-%m-%d %H:%M:%S"),
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
            "endTime": datetime.now(pytz.timezone("Europe/London"))
            .replace(tzinfo=None)
            .strftime("%Y-%m-%d %H:%M:%S"),
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
    ]
