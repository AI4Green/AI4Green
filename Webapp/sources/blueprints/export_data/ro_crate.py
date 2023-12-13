# from datetime import datetime
# from itertools import zip_longest
# from typing import Any, Dict, Generator, List, Tuple, Union
#
# import pytz
# from flask import requests
# from sources import models
#
#
# def get_author_json(reaction: models.Reaction) -> Dict:
#     person = reaction.creator_person
#     if person.user:
#         author = {
#             "@id": f"./author/{person.id}",
#             "@type": "Person",
#             "name": person.user.fullname,
#             "email": person.user.email,
#             "workgroup_role": person.user,
#         }
#     else:
#         author = {
#             "@id": f"./author/{person.id}",
#             "@type": "Person",
#             "name": "deleted",
#             "email": "deleted",
#             "workgroup_role": person.user,
#         }
#     return author
#
#
# def get_meta_data_json() -> Dict:
#     """
#     Returns the contents to make the metadata json file.
#     """
#     return {
#         "@id": "ro-crate-metadata.json",
#         "@type": "CreativeWork",
#         "about": {"@id": "./"},
#         "conformsTo": {"@id": "https://w3id.org/ro/crate/1.1"},
#         "dateCreated": datetime.now(pytz.timezone("Europe/London")).replace(
#             tzinfo=None
#         ),
#         "sdPublisher": {
#             "@type": "Organization",
#             "name": "AI4Green",
#             # "logo": "https://www.elabftw.net/img/elabftw-logo-only.svg",
#             # "slogan": "A free and open source electronic lab notebook.",
#             "url": "https://www.ai4green.app",
#             "parentOrganization": {
#                 "@type": "Organization",
#                 "name": "University of Nottingham",
#                 # "logo": "https://www.deltablot.com/img/logos/deltablot.svg",
#                 # "slogan": "Open Source software for research labs.",
#                 "url": "https://www.nottingham.ac.uk/",
#             },
#         },
#         "version": "1.0",
#     }
#
#
# def get_comment_id_list(reaction: models.Reaction) -> List[Dict]:
#     """
#     Returns a list of dictionaries, 1 per comment containing the id of the comment
#     """
#     comment_id_list = []
#     for comment in reaction.addenda:
#         # make an id from the work
#         comment_id_list.append(
#             {"@id": f"comment-id/{reaction.reaction_id}/{comment.id}"}
#         )
#     return comment_id_list
#
#
# def get_has_parts_id_list(reaction: models.Reaction) -> List[Dict]:
#     """
#     Returns a list of dictionaries, 1 per hasPart / file attachments and other files generated for export.
#     """
#     has_parts_list = []
#     for file in reaction.file_attachments:
#         has_parts_list.append({"@id": f"file-id/{reaction.reaction_id}/{file.id}"})
#     # include any MOL, SDF or RDF files. - id would follow same pattern but swap file.id for file extension.
#     return has_parts_list
#
#
# def get_file_details(file_id: str, reaction: models.Reaction) -> Tuple[Dict, str]:
#     """
#     Returns a string describing the file from its extension type
#     """
#     file_details = {}
#     name = ""
#     return file_details, name
#
#
# def get_has_parts_definitions(
#     reaction: models.Reaction, id_list: List[Dict]
# ) -> Generator[Dict[str]]:
#     """
#     Returns a Tuple of dictionaries containing the definitions of the hasParts which were identified earlier
#     """
#     definition_list = []
#     for file_id, file in zip_longest(id_list, reaction.file_attachments):
#         if file:
#             file_details = file.file_details
#             name = file.name
#         else:
#             name, file_details = get_file_details(file_id["@id"], reaction)
#
#         definition_list.append(
#             {
#                 "@id": file_id["@id"],
#                 "@type": "File",
#                 "description": "",  # Optional?
#                 "name": name,
#                 "encodingFormat": file_details["mimetype"],
#                 "contentSize": file_details["size"],
#                 "sha256": file_details[
#                     "sha256"
#                 ],  # TODO - need to make process to generate this somewhere.
#             }
#         )
#     return (definition for definition in definition_list)
#
#     #  also remember one of: [ .mol, .SDF, .RDF], csv with headings required to fill in the table, experimental.txt, pdf summary of rxn if not already
#
#
# def get_comment_definitions(reaction: models.Reaction, id_list):
#     comment_list = []
#     for comment_id, comment in zip(reaction, id_list):
#         comment_list.append(
#             {
#                 # comment properties
#             }
#         )
#     return (comment for comment in comment_list)
#
#
# def make_files(reaction):
#     """
#     Saves files to tmp storage
#     """
#     # make .rxn
#     # make experimental.txt
#     # make reaction.csv # data to reload reaction table
#
#
# def reaction_folder(workgroup_name: str, workbook_name: str, reaction: models.Reaction):
#     author_json = get_author_json(reaction)
#     make_files(reaction)
#
#     {
#         "@id": workgroup_name + "/" + workbook_name + "/" + reaction.reaction_id,
#         "@type": "Dataset",
#         "name": reaction.name,
#         "identifier": "20230928-ec8566c655af7310cd6c1ff77882d6a916191ead",
#         "author": {"@id": f"./author/{reaction.creator_person.id}"},
#         "dateCreated": reaction.time_of_creation,
#         "dateModified": reaction.time_of_update,
#         "comment": get_comment_id_list(reaction),
#         "hasPart": get_has_parts_id_list(reaction),  # files
#         "url": f"{requests.base_url}/{reaction.workgroup}/{reaction.workbook}/{reaction.reaction_id}/no",  # may remove no in future
#     }
#
#     # file json example
#
#     # {
#     #   "@id": "./2025-03-03 - Testing-the-eLabFTW-lab-notebook - b6db898e/example.png",
#     #   "@type": "File",
#     #   "description": "",
#     #   "name": "example.png",
#     #   "alternateName": "6b/6befead611faa1c99829ee811f5973e308f28cc8faea98594efc9ecf637bb62db3da3172e5c53e23f03bec7750b7ead5a32542b099544f8ee2902e4d66a99583.png",
#     #   "contentSize": 40959,
#     #   "sha256": "f2780bb5a8883f8c89a6da1c5bf4604f7cb44bdcb093e16484eb81c1bcb1f580"
#     # },
#     #
#     #
#     #
#     #
#     #
#     # {
#     #     "@id": "./2025-03-03 - Testing-the-eLabFTW-lab-notebook - b6db898e",
#     #     "@type": "Dataset",
#     #     "author": {
#     #         "@id": "person://a0568e4aa35eec2cf3962db57eec039a205d3b7d7c7132388af976dfaf8f2694?hash_algo=sha256"
#     #     },
#     #     "dateCreated": "2023-09-28T03:44:54+02:00",
#     #     "dateModified": "2023-10-02T14:51:49+02:00",
#     #     "identifier": "20230928-b6db898e231dbc2234227d19656d7e19a77654b7",
#     #     "comment": [
#     #         {
#     #             "@id": "comment://2023-09-28T03%3A44%3A54%2B02%3A00"
#     #         }
#     #     ],
#     #     "keywords": [
#     #         "demo",
#     #         "test"
#     #     ],
#     #     "name": "Testing the eLabFTW lab notebook",
#     #     "text": "<h1>Goal</h1>\n<p>Test the software.</p>\n<h1>Procedure</h1>\n<p>Click everywhere and explore everything.</p>\n<h1>Results</h1>\n<p>It's really nice, I think I'll adopt it for our lab.</p>",
#     #     "url": "https://elab.local:3148/experiments.php?mode=view&id=264",
#     #     "hasPart": [
#     #         {
#     #             "@id": "./2025-03-03 - Testing-the-eLabFTW-lab-notebook - b6db898e/export-elabftw.json"
#     #         },
#     #         {
#     #             "@id": "./2025-03-03 - Testing-the-eLabFTW-lab-notebook - b6db898e/example.png"
#     #         },
#     #         {
#     #             "@id": "./2025-03-03 - Testing-the-eLabFTW-lab-notebook - b6db898e/example.pdf"
#     #         }
#     #     ],
#     #     "mentions": [
#     #         {
#     #             "@id": "https://elab.local:3148/database.php?mode=view&id=51"
#     #         },
#     #         {
#     #             "@id": "https://elab.local:3148/database.php?mode=view&id=50"
#     #         },
#     #         {
#     #             "@id": "https://elab.local:3148/experiments.php?mode=view&id=256"
#     #         },
#     #         {
#     #             "@id": "https://elab.local:3148/experiments.php?mode=view&id=257"
#     #         },
#     #         {
#     #             "@id": "https://elab.local:3148/experiments.php?mode=view&id=258"
#     #         }
#     #     ]
#     # },
#
#     #
#     #
#     #
#     #
#     #
#     #
#     #
#     #
#     #
#     #
#     #
#     #
#     # {
#     #     "@id": "./some-unique-id/23",
#     #     "@type": "Dataset",
#     #     "author": {
#     #         "@id": "./some-author-id/44"
#     #     },
#     #     "dateCreated": "2023-09-23T01:02:26+02:00",
#     #     "dateModified": "2023-09-27T23:02:44+02:00",
#     #     "comment": [
#     #         {
#     #             "@id": "./some-comment-id/91"
#     #         }
#     #     ],
#     # },
#     # {
#     #     "@id": "./some-comment-id/91",
#     #     "@type": "Comment",
#     #     "dateCreated": "2023-09-23T01:02:26+02:00",
#     #     "text": "This is the content of the comment.",
#     #     "author": {
#     #         "@id": "./some-author-id/44"
#     #     }
#     # },
#     # {
#     #     "@id": "./some-author-id/44"
#     #            "@type": "Person",
#     # "familyName": "Tapie",
#     # "givenName": "Bernard"
#     # }
#     # }
#     #
#     #
