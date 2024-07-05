import base64
import io
import json
import os.path
import uuid
from datetime import datetime
from typing import Dict, List
from urllib.parse import quote

import pytz
import sources.services.data_export.export
from bs4 import BeautifulSoup
from flask import abort, url_for
from sources import models, services


class ELNFileExport:
    """
    Class for making an ELN File export consisting of a directory per reaction and a metadata json at the root
    describing the contents: https://github.com/TheELNConsortium/TheELNFileFormat
    """

    def __init__(self, data_export_request: models.DataExportRequest):
        """
        Creates an instance of the ELNFileExport class and names the root variable
        Args:
            data_export_request - the data export request we are making the eln file export for.
        """
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
        """The main function that calls the required functions to make the ELN file"""
        self._define_components()
        self._make_ro_crate_metadata_json()
        self._make_reaction_folders()
        self._remove_temp_blobs()

    def _remove_temp_blobs(self):
        """Deletes the temporary blobs used to create the export file"""
        container_client = (
            services.file_attachments.create_or_get_existing_container_client(
                self.blob_service_client, self.container_name
            )
        )
        blob_list = container_client.list_blobs(name_starts_with="temp/")

        # delete each blob in the temporary container
        for blob in blob_list:
            blob_client = self.blob_service_client.get_blob_client(
                container=self.container_name, blob=blob["name"]
            )
            blob_client.delete_blob()
            if blob_client.exists():
                abort(401)

    def _make_reaction_folders(self):
        """Makes a subdirectory in the .ELN zip export for each reaction being exported"""
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
        """Defines the components being described in the metadata json"""
        self._define_authors()
        for reaction in self.data_export_request.reactions:
            export_rxn = ELNExportReaction(reaction, self)
            export_rxn.make_export_files()
            export_rxn.get_files()
            export_rxn.get_parts()
            self.reaction_list.append(export_rxn)
            # get the root parts
            self.root_parts.append(export_rxn.defined_dataset)

    def _make_ro_crate_metadata_json(self):
        """Makes the metadata json file that describes the export contents."""
        ro_crate_graph = []
        ro_crate_graph += get_constant_ro_crate_start()
        # update the ro-crate with export specific data
        ro_crate_graph.append(self._describe_root_directory())

        for rxn in self.reaction_list:
            ro_crate_graph.append(rxn.defined_dataset)

        for rxn in self.reaction_list:
            [ro_crate_graph.append(file) for file in rxn.defined_files]

        for rxn in self.reaction_list:
            [ro_crate_graph.append(comment) for comment in rxn.defined_comments]

        ro_crate_graph += [author for author in self.defined_authors]

        # define the contents of the metadata json
        ro_crate_metadata_contents = {
            "@context": "https://w3id.org/ro/crate/1.1/context",
            "@graph":
            # the graph is a list of dictionaries defining the elements and their structure
            ro_crate_graph,
        }
        # save as a file to memory then upload the blob to azure.
        self.memory_file = io.StringIO()

        with self.memory_file as f:
            json.dump(ro_crate_metadata_contents, f)
            self.file_contents = bytearray(f.getvalue(), "utf-8")
            self.content_size = f.seek(0, os.SEEK_END)

        sources.services.data_export.export.save_blob(
            self.container_name,
            self.root + "/" + "ro-crate-metadata.json",
            self.file_contents,
        )

    def _describe_root_directory(self):
        """Describes the root directory by listing the reaction IDs in hasPart for the metadata json"""
        return {
            "@id": "./",
            "@type": "Dataset",
            "hasPart": [{"@id": part["@id"]} for part in self.root_parts],
        }

    def _define_authors(self):
        """List of all authors of the exported reactions"""
        [
            self.defined_authors.append(self._get_author(reaction))
            for reaction in self.data_export_request.reactions
        ]
        # remove duplicates
        self.defined_authors = [
            i
            for n, i in enumerate(self.defined_authors)
            if i not in self.defined_authors[:n]
        ]

    @staticmethod
    def _get_author(reaction: models.Reaction) -> Dict:
        """
        Gets the author data from the reaction creator property
        Args:
            reaction: the reaction being exported
        Returns:
            Details of the reaction author in the format required for the metadata JSON
        """
        person = reaction.creator_person
        if person.user:
            author = {
                "@id": f"./author/{person.id}",
                "@type": "Person",
                "fullName": person.user.fullname,
                "emailAddress": person.user.email,
            }
        # if the user has deleted their account then no personal details are included, only a unique ID is provided
        else:
            author = {
                "@id": f"./author/{person.id}",
                "@type": "Person",
                "fullName": "deleted",
                "emaiLAddress": "deleted",
            }
        return author


class ELNExportReaction:
    """Class for getting the reaction data required for an ELN File export."""

    def __init__(self, reaction: models.Reaction, eln_export: ELNFileExport):
        """
        Creates a new instance of the ELNExportReaction class
        Args:
            reaction - the reaction being exported
            eln_export - the eln file export that the exported reaction is a part of
        """

        self.reaction = reaction
        self.container_name = eln_export.data_export_request.uuid
        self.root = eln_export.root
        self.data_export_request = eln_export.data_export_request
        self.time = eln_export.time
        self.files = []
        self.defined_comments = []
        self.defined_files = []
        self.defined_dataset = []
        self.metadata = services.data_export.metadata.ReactionMetaData(reaction)
        self.rxn_metadata = self.metadata.get_dict()
        self.summary_soup = None

    def make_export_files(self):
        """Makes files which are included with the ELN File export. Includes: .rdf"""
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
            "uuid": blob_name,
            "ro-crate-id": file_uuid,
        }
        self.files.append(rdf_dict)

    def get_files(self):
        """Gets the file details for the reaction being exported"""
        # get pre existing files
        pre_existing_files = [
            x.to_export_dict() for x in self.reaction.file_attachments
        ]
        for file in pre_existing_files:
            file["made_for_export"] = False
        self.files += pre_existing_files
        self._make_file_display_names_unique()

    def _make_file_display_names_unique(self):
        """Change display names to ensure they are all unique - required for the metadata json ID"""
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
        """One reaction is its own dataset within the ELN file export and is defined here."""
        self.defined_dataset = {
            "@id": f"./{self.reaction.reaction_id}",
            "@type": "Dataset",
            "name": self.reaction.name,  # experiment title
            "author": {"@id": f"./author/{self.reaction.creator_person.id}"},
            "dateCreated": self.reaction.time_of_creation.strftime("%Y-%m-%d %H:%M:%S"),
            "dateModified": self.reaction.time_of_update.strftime("%Y-%m-%d %H:%M:%S"),
            "comment": [comment["@id"] for comment in self.defined_comments],
            "hasPart": [{"@id": file["@id"]} for file in self.defined_files],
            "description": self.reaction.description,
            "text": self._get_summary_html(),  # todo restore after debugging
            "encodingFormat": "text/html",
            "url": url_for("main.index", _external=True)
            + quote(
                f"/{self.reaction.workbook.WorkGroup.name}/{self.reaction.workbook.name}/{self.reaction.reaction_id}/no"
            ),
            # may remove no in future
        }

    def _define_comments(self):
        """Gets all the comments associated with a reaction and defines them in a format for the metadata json"""
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
        """Gets all the file details associated with a reaction and defines them in a format for the metadata json"""
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

    @staticmethod
    def _get_file_description(file: Dict) -> str:
        """Gets the description for a file attachment of a reaction"""
        if file["made_for_export"] is True:
            export_file_type_description_dict = {
                "RDF": "A Chemistry specific format with serialized metadata and molecules represented in a "
                "Chemical Table format and assigned as reactants, products, or agents"
            }
            description = export_file_type_description_dict[file["filetype"]]
        else:
            description = "Reaction data uploaded by user"
        return description

    def _extract_input_values_to_td(self):
        """
        Replaces <td><input ... value="..."></td> patterns with <td>value</td> in the summary_soup.
        """
        # check all td tags (table cells) and move any input tag values into text within the td tag
        for td in self.summary_soup.find_all("td"):
            # Check if there is an <input> element inside the <td>
            if td.find("input"):
                input_tag = td.find("input")
                # Ensure the <input> tag has a 'value' attribute and is not a radio or checkbox
                if (
                    input_tag.has_attr("value")
                    and input_tag["value"]
                    and not input_tag.attrs.get("type") in ["radio", "checkbox"]
                ):
                    # Extract the value from the <input> field and then remove it from the tree
                    input_value = input_tag["value"]
                    # Create a new text node with the extracted value
                    new_string = BeautifulSoup(input_value, "html.parser").text

                    # Add the new text node to the <td> while preserving existing content and child elements
                    if td.contents:
                        # Iterate over existing content and add new text node
                        for content in td.contents:
                            # If it's a NavigableString or whitespace, add the new string
                            if content.name is None or content.name == "br":
                                content.replace_with(new_string)
                            else:
                                # Append new text node before or after elements as needed
                                td.insert(td.contents.index(content) + 1, new_string)
                    else:
                        # If there are no existing contents, just set the new text value
                        td.string = new_string

    def _remove_unselected_options(self, keep_dict: Dict[str, str]):
        """
        Removes unselected <option> elements from <select> tags, keeping only the specified options.

        Args:
            keep_dict (Dict[str, str]): A dictionary mapping select ids to the text of the option to keep.
        """
        for select in self.summary_soup.find_all("select", {"id": True}):
            select_id = select["id"]
            if select_id in keep_dict:
                keep_text = keep_dict[select_id]
                for option in select.find_all("option"):
                    if option.text.strip() != keep_text:
                        option.decompose()

    def _remove_unchecked_checkboxes(self, id_list: List[str]) -> None:
        """
        Removes unchecked checkboxes from <td> tags based on a list of IDs.

        Args:
            id_list (List[str]): A list of checkbox ids to keep.
        """
        checkboxes = self.summary_soup.find_all("input", attrs={"type": "checkbox"})
        for checkbox in checkboxes:
            if checkbox.attrs["id"] not in id_list:
                checkbox.parent.decompose()

    def _update_values(self, value_dict: Dict[str, str]):
        for key, value in value_dict.items():
            # find either an input inside a td that matches the id or an input that matches the id itself.
            input_tag = (
                self.summary_soup.find("td", id=key).find("input", {"type": "number"})
                if self.summary_soup.find("td", id=key)
                else self.summary_soup.find("input", id=key)
            )
            if input_tag:
                input_tag["value"] = value  # Set the new value

    def _get_summary_html(self):
        """Creates a html summary of the pdf from the summary_to_print value and adapts for ELN import"""
        # get the summary table html and required data from the database object
        reaction_data = json.loads(self.reaction.reaction_table_data)
        summary_table_data = json.loads(self.reaction.summary_table_data)
        summary_table_html = summary_table_data["summary_to_print"]
        self.summary_soup = BeautifulSoup(summary_table_html, "html.parser")

        # make a dictionary with html element ids and their corresponding values to update the html with
        main_product_idx = int(reaction_data["main_product"]) - 1
        theoretical_yield = reaction_data["product_masses"][main_product_idx]
        id_value_dict = {
            "js-ae-cell": self.rxn_metadata["atom_economy"],
            "js-me-cell": self.rxn_metadata["mass_efficiency"],
            "js-yield-cell": self.rxn_metadata["percent_yield"],
            "js-percentage-yield": self.rxn_metadata["percent_yield"],
            "js-conversion-cell": self.rxn_metadata["conversion"],
            "js-selectivity-cell": self.rxn_metadata["selectivity"],
            "js-unreacted-reactant-mass": summary_table_data["unreacted_reactant_mass"],
            "js-real-product-mass": summary_table_data["real_product_mass"],
            "js-percentage-yield-cell": self.rxn_metadata["percent_yield"],
            "js-product-rounded-mass"
            + reaction_data["main_product"]: theoretical_yield,
            "js-temperature-cell": self.rxn_metadata["temperature"],
            "js-risk-score": self.metadata.exported_from_js["risk-score"],
            "js-researcher": "Researcher:" + summary_table_data["researcher"],
            "js-supervisor": "Supervisor:" + summary_table_data["supervisor"],
        }
        self._update_values(id_value_dict)

        # make a dictionary with the dropdowns and their selected values then remove the non-selected options
        dropdowns_value_dict = {
            "js-elements": self.rxn_metadata["element_sustainability"],
            "js-batch-flow": self.rxn_metadata["batch_or_flow"],
            "js-isolation": self.rxn_metadata["isolation_method"],
            "js-catalyst": self.rxn_metadata["catalyst_use"],
            "js-recovery": self.rxn_metadata["catalyst_recovery"],
        }
        self._remove_unselected_options(dropdowns_value_dict)

        # remove unselected checkboxes from the html
        checked_boxes = summary_table_data["radio_buttons"]
        self._remove_unchecked_checkboxes(checked_boxes)

        # update the href from a relative path to an absolute url
        a_tag = self.summary_soup.find("a", href="/element_sustainability")
        if a_tag:
            a_tag["href"] = "https://ai4green.app/element_sustainability"

        # get a list of checked risk boxes and make these bold in the html
        risk_boxes = [
            "slight",
            "serious",
            "major",
            "lowLikelihood",
            "possible",
            "frequentOccur",
            "individual",
            "localLabs",
            "buildingWide",
        ]
        checked_risk_boxes = [
            risk_box for risk_box in risk_boxes if risk_box in checked_boxes
        ]
        self._make_checked_risk_boxes_bold(checked_risk_boxes)

        # finally, update the html to replace <td><input value='x'></td> with <td>x</td>
        self._extract_input_values_to_td()
        summary_table_html = str(self.summary_soup)
        summary_html = f"<h2>{self.reaction.reaction_id}</h2><div>Reaction Smiles: {self.reaction.reaction_smiles}{summary_table_html}</div>"
        return summary_html

    def _make_checked_risk_boxes_bold(self, checked_risk_boxes: List[str]):
        """
        Makes the text of labels bold for the checkboxes specified in the checked_risk_boxes list.

        This method finds all <label> tags in the summary_soup and checks if the 'for' attribute matches any id in the checked_risk_boxes list.
        If there is a match, it wraps the label text in <b> tags to make it bold.

        Args:
            checked_risk_boxes (List[str]): A list of checkbox ids for which the associated label text should be made bold.
        """
        # Find all label tags
        for box in checked_risk_boxes:
            label = self.summary_soup.find("label", attrs={"for": box})
            # Check if the label was found
            if label:
                # Make the text bold using BeautifulSoup
                bold_tag = self.summary_soup.new_tag("b")
                bold_tag.string = label.text
                label.clear()
                label.append(bold_tag)


def get_constant_ro_crate_start() -> List[Dict]:
    """Returns data which is consistent across all ELN file exports"""
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
    """
    Gets the root directory name for the eln file export.
    Args:
        data_export_request - the request being exported
    Returns:
        the root name as a string. Combines workbook name and current date.
    """
    return (
        data_export_request.workbooks[0].name
        + "-"
        + data_export_request.time_of_request.strftime("%Y-%m-%d")
        + "-export"
    )
