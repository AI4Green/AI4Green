# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
# This module receives data from the reaction table
# via POST request and renders the summary template
# """
# import ast
# import json
#
# import rdkit
# from flask import Response, abort, jsonify, render_template, request
# from flask_login import login_required
# from sources import auxiliary, db, models, services
#
# from . import summary_bp
#
#
# # Processing data from the reaction table and creating summary with H&S report
# @summary_bp.route("/_summary", methods=["POST"])
# def summary() -> Response:
#     # must be logged in
#     """This function receives the reaction information from browser, calculates
#     green metrics, gives hazard information, and renders the summary"""
#     if not (
#         str(request.form["demo"]) == "demo" or str(request.form["tutorial"]) == "yes"
#     ):
#         # check user permission
#         workgroup_name = str(request.form["workgroup"])
#         workbook_name = str(request.form["workbook"])
#         if not auxiliary.security_member_workgroup_workbook(
#             workgroup_name, workbook_name
#         ):
#             abort(401)
#     # check there is data to get and if so get it
#     summary_table_data = str(request.form["js_summary_table_data"])
#     if summary_table_data != "no data":
#         summary_table_data = ast.literal_eval(summary_table_data)
#
#     category_rate = {
#         "": 0,
#         "L": 1,
#         "M": 2,
#         "H": 3,
#         "VH": 4,
#     }
#
#     # reactant_mass_factor = mass_factor[mass_unit]
#     # product_mass_factor = mass_factor[product_mass_unit]
#
#     unit_data = services.summary.get_unit_data()
#     reactant_data = services.summary.get_reactant_data()
#     reagent_data = services.summary.get_reagent_data()
#     solvent_data = services.summary.get_solvent_data()
#     product_data = services.summary.get_product_data()
#
#     # check all the requirement information has been typed into the reaction table
#     check_results = services.summary.check_required_data_is_present(
#         reactant_data, reagent_data, solvent_data, product_data
#     )
#     if check_results != "checks_successful":
#         return jsonify({"summary": check_results})
#
#     sustainability_data = services.summary.get_sustainability_metrics()
#     # FROM HERE TODO
#
#     ae = (
#         round(
#             100
#             * float(product_molecular_weights[main_prod_idx])
#             / (reactant_molecular_weight_sum + reagent_molecular_weight_sum),
#             1,
#         )
#         if (reactant_molecular_weight_sum + reagent_molecular_weight_sum) > 0
#         else 0
#     )  # atom economy
#     ae_flag = auxiliary.metric_flag(ae)  # atom economy color flag
#     # get reaction smiles and convert to a list of smiles involved in reaction
#     reaction_smiles = str(request.form["reactionSmiles"]).split(" |")[0]
#     reaction_smiles_ls = reaction_smiles.replace(">>", ".").split(".")
#
#     # get reagent and solvent smiles
#     reagent_smiles_ls = services.all_compounds.get_smiles_list(reagent_primary_keys_ls)
#     solvent_smiles_ls = services.all_compounds.get_smiles_list(solvent_primary_keys_ls)
#
#     # reaction smiles already has reactant and product smiles
#     full_reaction_smiles_ls = [
#         x for x in reaction_smiles_ls + reagent_smiles_ls + solvent_smiles_ls if x
#     ]
#     full_reaction_smiles_ls = [x for x in full_reaction_smiles_ls if x]
#     # get list of elements in reaction from the smiles list by converting to atoms and then chemical symbols
#     element_symbols = set()
#     for component in full_reaction_smiles_ls:
#         mol = rdkit.Chem.MolFromSmiles(component)
#         for atom in mol.GetAtoms():
#             symbol = atom.GetSymbol()
#             element_symbols.add(symbol)
#     # get sustainability colour and set element_sustainability variable according to least sustainable element present
#     element_sustainability_set = set(
#         y.colour
#         for y in [
#             db.session.query(models.Element.colour)
#             .filter(models.Element.symbol == symbol)
#             .first()
#             for symbol in element_symbols
#         ]
#     )
#     element_sustainability = "-select-"
#     element_sustainability_flag = "hazard-reset-hazard"
#     if "red" in element_sustainability_set:
#         element_sustainability = "5-50 years"
#         element_sustainability_flag = "hazard-hazardous"
#     elif "yellow" in element_sustainability_set:
#         element_sustainability = "50-500 years"
#         element_sustainability_flag = "hazard-warning"
#     elif "lime" in element_sustainability_set:
#         element_sustainability = "+500 years"
#         element_sustainability_flag = "hazard-acceptable"
#
#     # To here
#
#     # Hazard summary
#     # Reactant hazards
#     (
#         reactant_most_severe_hazard_numerical_rating,
#         reactant_hazard_sentences,
#         reactant_hazard_ratings,
#         reactant_hazard_colors,
#         reactant_exposure_potentials,
#         reactant_risk_ratings,
#         reactant_risk_colors,
#     ) = services.hazard_code.get_multiple_compounds_data(
#         reagent_data["reactant_hazards"], reactant_data["reactant_physical_forms"]
#     )
#
#     # Reagent hazards
#     (
#         reagent_most_severe_hazard_numerical_rating,
#         reagent_hazard_sentences,
#         reagent_hazard_ratings,
#         reagent_hazard_colors,
#         reagent_exposure_potentials,
#         reagent_risk_ratings,
#         reagent_risk_colors,
#     ) = services.hazard_code.get_multiple_compounds_data(
#         reagent_data["reagent_hazards"], reagent_data["reagent_physical_forms"]
#     )
#     # solvent hazards
#     (
#         solvent_most_severe_hazard_numerical_rating,
#         solvent_hazard_sentences,
#         solvent_hazard_ratings,
#         solvent_hazard_colors,
#         solvent_exposure_potentials,
#         solvent_risk_ratings,
#         solvent_risk_colors,
#     ) = services.hazard_code.get_multiple_compounds_data(
#         solvent_data["solvent_hazards"], solvent_data["solvent_physical_forms"]
#     )
#
#     # Product hazard
#     (
#         product_most_severe_hazard_numerical_rating,
#         product_hazard_sentences,
#         product_hazard_ratings,
#         product_hazard_colors,
#         product_exposure_potentials,
#         product_risk_ratings,
#         product_risk_colors,
#     ) = services.hazard_code.get_multiple_compounds_data(
#         product_data["product_hazards"], product_data["product_physical_forms"]
#     )
#
#     most_severe_hazard_numerical_rating = (
#         reactant_most_severe_hazard_numerical_rating
#         + reagent_most_severe_hazard_numerical_rating
#         + solvent_most_severe_hazard_numerical_rating
#         + product_most_severe_hazard_numerical_rating
#     )
#     max_most_severe_hazard_numerical_rating = int(
#         max(most_severe_hazard_numerical_rating)
#     )  # max total hazard rate
#     risk_rating = list(category_rate.keys())[
#         max_most_severe_hazard_numerical_rating
#     ]  # resulting hazard rating
#     risk_color = (
#         "hazard-hazardous" if risk_rating == "VH" else "hazard-reset-hazard"
#     )  # color code for the hazard rating
#
#     # Solvent flags
#     flag_rate = {
#         1: "hazard-highly-hazardous",
#         2: "hazard-hazardous",
#         3: "hazard-warning",
#         4: "hazard-acceptable",
#         5: "non-chem21",
#     }  # flag rate dictionary
#     solvent_flags = []  # solvent flag list
#     if solvents[0]:
#         for solvent in solvents:
#             solvent_flag = (
#                 db.session.query(models.Solvent.flag)
#                 .filter(models.Solvent.name == solvent)
#                 .first()
#             )
#             solvent_flag = solvent_flag[0] if solvent_flag else None
#             try:
#                 solvent_flags.append(
#                     flag_rate[solvent_flag]
#                 )  # appends solvent flag to their list
#             except KeyError:
#                 solvent_flags.append(5)
#
#     if number_of_solvents == "0":  # if no solvents have been chosen
#         number_of_solvents = 1  # then it shows only one empty cell
#         solvents = [" "]
#
#     services.summary.make_rxn_file()
#
#     # if product mass and reactant mass sum are calculated, then it forms a summary table
#     if product_data and reactant_data:
#         summary_table = render_template(
#             "_summary_table.html",
#             amount_unit=unit_data["amount_unit"],
#             volume_unit=unit_data["volume_unit"],
#             mass_unit=unit_data["mass_unit"],
#             solvent_volume_unit=unit_data["solvent_volume_unit"],
#             product_mass_unit=unit_data["product_mass_unit"],
#             reactants=reactant_data["reactants"],
#             reactant_primary_keys=reactant_data["reactant_primary_keys"],
#             reagent_primary_keys=reagent_data["reagent_primary_keys_str"],
#             reagents=reagent_data["reagents"],
#             reagent_table_numbers=reagent_data["reagent_table_numbers"],
#             reagent_molecular_weights=reagent_data["reagent_molecular_weights"],
#             reagent_densities=reagent_data["reagent_densities"],
#             reagent_concentrations=reagent_data["reagent_concentrations"],
#             reagent_equivalents=reagent_data["reagent_equivalents"],
#             reagent_hazards=reagent_data["reagent_hazards"],
#             reagent_amounts=reagent_data["reagent_amounts"],
#             rounded_reagent_amounts=reagent_data["rounded_reagent_amounts"],
#             reagent_volumes=reagent_data["reagent_volumes"],
#             rounded_reagent_volumes=reagent_data["rounded_reagent_volumes"],
#             reagent_masses=reagent_data["reagent_masses"],
#             rounded_reagent_masses=reagent_data["rounded_reagent_masses"],
#             solvents=solvent_data["solvents"],
#             solvent_volumes=solvent_data["solvent_volumes"],
#             solvent_table_numbers=solvent_data["solvent_table_numbers"],
#             solvent_flags=solvent_data["solvent_flags"],
#             products=product_data["products"],
#             product_table_numbers=product_data["product_table_numbers"],
#             reactant_molecular_weights=reactant_data["reactant_molecular_weights"],
#             reactant_densities=reactant_data["reactant_densities"],
#             reactant_concentrations=reactant_data["reactant_concentrations"],
#             reactant_equivalents=reactant_data["reactant_equivalents"],
#             reactant_amounts=reactant_data["reactant_amounts"],
#             rounded_reactant_amounts=reactant_data["rounded_reactant_amounts"],
#             reactant_volumes=reactant_data["reactant_volumes"],
#             rounded_reactant_volumes=reactant_data["rounded_reactant_volumes"],
#             reactant_masses=reactant_data["reactant_masses"],
#             rounded_reactant_masses=reactant_data["rounded_reactant_masses"],
#             product_primary_keys=product_data["product_primary_keys"],
#             main_product_table_number=product_data["main_product_table_number"],
#             main_product_index=product_data["main_product_index"],
#             product_molecular_weights=product_data["product_molecular_weights"],
#             product_masses=product_data["product_masses"],
#             rounded_product_masses=product_data["rounded_product_masses"],
#             ae=product_data["ae"],
#             ae_flag=product_data["ae_flag"],
#             element_sustainability=product_data["element_sustainability"],
#             element_sustainability_flag=product_data["element_sustainability_flag"],
#             reactant_hazard_sentences=reactant_data["reactant_hazard_sentences"],
#             reactant_hazard_ratings=reactant_data["reactant_hazard_ratings"],
#             reactant_hazard_colors=reactant_data["reactant_hazard_colors"],
#             reactant_risk_colors=reactant_data["reactant_risk_colors"],
#             reactant_exposure_potentials=reactant_data["reactant_exposure_potentials"],
#             reactant_risk_ratings=reactant_data["reactant_risk_ratings"],
#             reagent_hazard_sentences=reagent_data["reagent_hazard_sentences"],
#             reagent_hazard_ratings=reagent_data["reagent_hazard_ratings"],
#             reagent_hazard_colors=reagent_data["reagent_hazard_colors"],
#             reagent_risk_colors=reagent_data["reagent_risk_colors"],
#             reagent_exposure_potentials=reagent_data["reagent_exposure_potentials"],
#             reagent_risk_ratings=reagent_data["reagent_risk_ratings"],
#             solvent_primary_keys=solvent_data["solvent_primary_keys_str"],
#             solvent_hazard_sentences=solvent_hazard_sentences,
#             solvent_hazard_ratings=solvent_hazard_ratings,
#             solvent_exposure_potentials=solvent_exposure_potentials,
#             solvent_risk_ratings=solvent_risk_ratings,
#             solvent_hazard_colors=solvent_hazard_colors,
#             solvent_risk_colors=solvent_risk_colors,
#             product_hazard_sentences=product_hazard_sentences,
#             product_hazard_ratings=product_hazard_ratings,
#             product_exposure_potentials=product_exposure_potentials,
#             product_risk_ratings=product_risk_ratings,
#             product_hazard_colors=product_hazard_colors,
#             product_risk_colors=product_risk_colors,
#             risk_rating=risk_rating,
#             risk_color=risk_color,
#             number_of_solvents=number_of_solvents,
#             summary_table_data=json.dumps(summary_table_data),
#         )
#         return jsonify({"summary": summary_table})
#     else:
#         pass
#     return jsonify(
#         {"summary": "Ensure you have entered all the necessary information!"}
#     )  # otherwise it shows this message
#
#
# @summary_bp.route("/element_sustainability", methods=["POST", "GET"])
# @login_required
# def element_sustainability() -> Response:
#     # must be logged in
#     return render_template("element_sustainability.html")
#
#
# @summary_bp.route("/pdf", methods=["POST", "GET"])
# @login_required
# def pdf():
#     """
#     Saves the autogenerated PDF for a reaction, overwriting the old one if present.
#     """
#     services.auth.reaction_files(permission_level="edit")
#     new_upload = services.file_attachments.UploadExperimentDataFiles(
#         request, autogenerated_file=True
#     )
#     new_upload.validate_files()
#     new_upload.remove_duplicate_autogenerated_summaries()
#     new_upload.remove_duplicate_autogenerated_summaries()
#     new_upload.save_validated_files()
#     return "", 204
#
#
# @summary_bp.route("/get_file_attachment_list", methods=["POST"])
# @login_required
# def get_file_attachment_list():
#     """
#     Gets a list of file attachments for the current reaction identified from the request data
#     """
#     services.auth.reaction_files(permission_level="view")
#     reaction = services.reaction.get_current_from_request()
#     file_attachments = sorted(
#         reaction.file_attachments, key=lambda x: not x.autogenerated
#     )
#     file_attachments_dict_list = [
#         {"name": x.display_name, "uuid": x.uuid, "autogenerated": x.autogenerated}
#         for x in file_attachments
#     ]
#     return jsonify({"file_attachments": file_attachments_dict_list})
