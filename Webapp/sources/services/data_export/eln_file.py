import io
import json
import os.path
import uuid
from datetime import datetime
from itertools import zip_longest
from typing import Dict, Generator, List, Optional, Tuple
from urllib.parse import quote

import pytz
import sources.services.data_export.export
from flask import abort, current_app, request, url_for
from rdkit.Chem import AllChem
from sources import models, services

# todo need to use './' in the json and 'root/' in the blob names
# solution use: ro-crate-id and blob-id


class ELNFileExport:
    # absolute basic example for now just include an rdf for a reaction

    def __init__(self, data_export_request: models.DataExportRequest):
        self.data_export_request = data_export_request
        self.container_name = data_export_request.uuid
        self.blob_service_client = (
            services.file_attachments.connect_to_azure_blob_service()
        )
        self.root = get_root_name(self.data_export_request)

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
        self._define_components()
        self._make_ro_crate_metadata_json()
        self._make_reaction_folders()
        self._remove_temp()

    def _remove_temp(self):
        container_client = (
            services.file_attachments.create_or_get_existing_container_client(
                self.blob_service_client, self.container_name
            )
        )
        blob_list = container_client.list_blobs(name_starts_with="temp/")

        # blob_list =
        for blob in blob_list:
            print(blob["name"])
            blob_client = self.blob_service_client.get_blob_client(
                container=self.container_name, blob=blob["name"]
            )
            blob_client.delete_blob()
            if blob_client.exists():
                abort(401)

    def _make_reaction_folders(self):
        for reaction in self.reaction_list:
            # for each reaction make a folder named after the dataset (workbook) with all belonging files included
            subdirectory = reaction.reaction.reaction_id + "/"
            # for each file copy into the new directory
            for file in reaction.files:
                # find the original blob with the container and the uuid
                source_blob_client = self.blob_service_client.get_blob_client(
                    container=file["container_name"], blob=file["uuid"]
                )
                # make a new blob in the temporary export container in the subdir with the display name
                new_blob = self.blob_service_client.get_blob_client(
                    container=self.container_name,
                    blob=self.root + "/" + subdirectory + file["display_name"],
                )
                new_blob.start_copy_from_url(source_blob_client.url)

    def _define_components(self):
        self._define_authors()
        for reaction in self.data_export_request.reactions:
            export_rxn = ELNExportReaction(reaction, self)
            export_rxn.get_files()
            export_rxn.get_parts()
            self.reaction_list.append(export_rxn)
            # get the root parts
            self.root_parts = export_rxn.defined_datasets

    def _make_ro_crate_metadata_json(self):
        ro_crate_graph = []
        # this part is common to all .eln files from AI4Green
        ro_crate_graph += get_constant_ro_crate_start()
        ro_crate_graph.append(self._describe_root_directory())

        ro_crate_graph += [
            dataset for rxn in self.reaction_list for dataset in rxn.defined_datasets
        ]

        ro_crate_graph += [
            file for rxn in self.reaction_list for file in rxn.defined_files
        ]

        ro_crate_graph += [
            comment for rxn in self.reaction_list for comment in rxn.defined_comments
        ]

        ro_crate_graph += [author for author in self.defined_authors]

        ro_crate_metadata_contents = {
            "@context": "https://w3id.org/ro/crate/1.1/context",
            "@graph":
            # the graph is a list of dictionaries defining the elements and their structure
            ro_crate_graph,
        }
        self.memory_file = io.StringIO()

        with self.memory_file as f:
            json.dump(ro_crate_metadata_contents, f)
            self.file_contents = bytearray(f.getvalue(), "utf-8")
            self.content_size = f.seek(0, os.SEEK_END)

        # save this json to azure in root dir
        # blob_client = self.blob_service_client.get_blob_client(self.container_name, "ro_crate_metadata.json")
        # blob_client.upload_blob()
        sources.services.data_export.export.save_blob(
            self.container_name,
            self.root + "/" + "ro-crate-metadata.json",
            self.file_contents,
        )

        # # testing the json looks right
        # with open("ro-crate-metadata.json", "w") as f:
        #     f.write(json.dumps(ro_crate_metadata_contents))

    def _describe_root_directory(self):
        return {
            "@id": "./",
            "@type": "Dataset",
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
        self.container_name = eln_export.data_export_request.uuid
        self.root = eln_export.root
        self.data_export_request = eln_export.data_export_request
        self.time = eln_export.time
        self.files = []
        self.defined_comments = []
        self.defined_files = []
        self.defined_datasets = []

    def _make_export_files(self):
        # for now we JUST make rdf
        file_uuid = str(uuid.uuid4())
        blob_name = "temp/" + file_uuid
        rdf = services.data_export.serialisation.ReactionDataFileExport(
            self.reaction, blob_name, self.container_name
        )
        rdf.save(extension=False)
        rdf_dict = {
            "filetype": "RDF",
            "made_for_export": True,
            "display_name": self.reaction.reaction_id + ".rdf",
            "container_name": rdf.container_name,
            "content_size": rdf.content_size,
            "sha256": rdf.file_hash,
            "mimetype": rdf.mime_type,
            "time_of_creation": self.time,
            # "temp_uuid": blob_name,
            "uuid": blob_name,
            "ro-crate-id": file_uuid
            # "blob_name": blob_name
        }
        self.files.append(rdf_dict)

    def get_files(self):
        # self._make_export_files()
        # get pre existing files
        pre_existing_files = [
            x.to_export_dict() for x in self.reaction.file_attachments
        ]
        for file in pre_existing_files:
            file["made_for_export"] = False
        self.files += pre_existing_files
        self._make_file_display_names_unique()

    def _make_file_display_names_unique(self):
        """Change display names to ensure they are all unique"""
        file_names = []
        for file in self.files:
            original_name, ext = os.path.splitext(file["display_name"])
            new_name = file["display_name"]
            count = 0
            while new_name in file_names:
                count += 1
                new_name = f"{original_name}_{count}{ext}"
            file["display_name"] = new_name
            file_names.append(new_name)

    def get_parts(self):
        """Define all parts for a reaction including the comments, files, and dataset"""
        self._define_comments()
        self._define_files()
        self._define_reaction_dataset()

    def _define_reaction_dataset(self):
        # each reaction has 1 part for the dataset (reaction) and n number of files for the reaction
        self.defined_datasets.append(
            {
                "@id": f"./{self.reaction.reaction_id}",
                "@type": "Dataset",
                "name": self.reaction.name,  # experiment title
                # "identifier": self.reaction.uuid, reaction has no uuid yet.
                "author": {"@id": f"./author/{self.reaction.creator_person.id}"},
                "dateCreated": self.reaction.time_of_creation.strftime(
                    "%Y-%m-%d %H:%M:%S"
                ),
                "dateModified": self.reaction.time_of_update.strftime(
                    "%Y-%m-%d %H:%M:%S"
                ),
                "comment": [comment["@id"] for comment in self.defined_comments],
                "hasPart": [
                    {"@id": file["@id"]} for file in self.defined_files
                ],  # files
                "description": "<h1>Chemistry Time</h1>",
                "text": "more chemistry",  # todo format as html for summary table + write up
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
                    "@id": f"./{self.reaction.reaction_id}/{file['display_name']}",
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
                "Chemical Table format and assigned as reactants, products, or agents"
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
                "name": "AI4Green",
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


def get_root_name(data_export_request: models.DataExportRequest) -> str:
    return (
        data_export_request.workbooks[0].name
        + "-"
        + data_export_request.time_of_request.strftime("%Y-%m-%d")
        + "-export"
    )
