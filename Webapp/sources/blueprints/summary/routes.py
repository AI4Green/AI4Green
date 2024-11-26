#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module receives data from the reaction table
via POST request and renders the summary template
"""
import ast
import json

from flask import Response, abort, jsonify, render_template, request
from flask_login import login_required
from sources import auxiliary, services

from . import summary_bp


# Processing data from the reaction table and creating summary with H&S report
@summary_bp.route("/_summary", methods=["POST"])
def summary() -> Response:
    # must be logged in
    """This function receives the reaction information from browser, calculates
    green metrics, gives hazard information, and renders the summary"""
    if not (
        str(request.form["demo"]) == "demo" or str(request.form["tutorial"]) == "yes"
    ):
        # check user permission
        workgroup_name = str(request.form["workgroup"])
        workbook_name = str(request.form["workbook"])
        if not auxiliary.security_member_workgroup_workbook(
            workgroup_name, workbook_name
        ):
            abort(401)

    unit_data = services.summary.get_unit_data(request.form)
    reactant_data = services.summary.get_reactant_data(request.form)
    reagent_data = services.summary.get_reagent_data(request.form)
    solvent_data = services.summary.get_solvent_data(request.form)
    product_data = services.summary.get_product_data(request.form)

    # check all the requirement information has been typed into the reaction table
    check_results = services.summary.check_required_data_is_present(
        reactant_data, reagent_data, solvent_data, product_data
    )
    # if the check fails we return a message to inform the user they are missing data
    if check_results != "checks successful":
        return jsonify({"summary": check_results})

    sustainability_data = services.sustainability.SustainabilityMetrics(
        reactant_data, reagent_data, solvent_data, product_data
    ).get_sustainability_metrics()

    risk_data = services.summary.get_risk_data(
        reactant_data, reagent_data, solvent_data, product_data
    )

    polymer_mode = request.form["polymerMode"]
    # change summary table in polymer mode
    if polymer_mode.lower() == "true":
        summary_table_html = "_polymer_summary_table.html"
    else:
        summary_table_html = "_summary_table.html"

    # if product mass and reactant mass sum are calculated, then it forms a summary table
    if product_data and reactant_data:
        summary_table = render_template(
            summary_table_html,
            amount_unit=unit_data["amount_unit"],
            volume_unit=unit_data["volume_unit"],
            mass_unit=unit_data["mass_unit"],
            solvent_volume_unit=unit_data["solvent_volume_unit"],
            product_mass_unit=unit_data["product_mass_unit"],
            reactants=reactant_data["reactants"],
            reactant_primary_keys=reactant_data["reactant_primary_keys_str"],
            reagent_primary_keys=reagent_data["reagent_primary_keys_str"],
            reagents=reagent_data["reagents"],
            reagent_table_numbers=reagent_data["reagent_table_numbers"],
            reagent_molecular_weights=reagent_data["reagent_molecular_weights"],
            reagent_densities=reagent_data["reagent_densities"],
            reagent_concentrations=reagent_data["reagent_concentrations"],
            reagent_equivalents=reagent_data["reagent_equivalents"],
            reagent_hazards=reagent_data["reagent_hazards"],
            reagent_amounts=reagent_data["reagent_amounts"],
            rounded_reagent_amounts=reagent_data["rounded_reagent_amounts"],
            reagent_volumes=reagent_data["reagent_volumes"],
            rounded_reagent_volumes=reagent_data["rounded_reagent_volumes"],
            reagent_masses=reagent_data["reagent_masses"],
            rounded_reagent_masses=reagent_data["rounded_reagent_masses"],
            solvents=solvent_data["solvents"],
            solvent_volumes=solvent_data["solvent_volumes"],
            solvent_table_numbers=solvent_data["solvent_table_numbers"],
            solvent_flags=sustainability_data["solvent_flags"],
            products=product_data["products"],
            product_table_numbers=product_data["product_table_numbers"],
            reactant_molecular_weights=reactant_data["reactant_molecular_weights"],
            reactant_densities=reactant_data["reactant_densities"],
            reactant_concentrations=reactant_data["reactant_concentrations"],
            reactant_equivalents=reactant_data["reactant_equivalents"],
            reactant_amounts=reactant_data["reactant_amounts"],
            rounded_reactant_amounts=reactant_data["rounded_reactant_amounts"],
            reactant_volumes=reactant_data["reactant_volumes"],
            rounded_reactant_volumes=reactant_data["rounded_reactant_volumes"],
            reactant_masses=reactant_data["reactant_masses"],
            rounded_reactant_masses=reactant_data["rounded_reactant_masses"],
            product_primary_keys=product_data["product_primary_keys_str"],
            main_product_table_number=product_data["main_product_table_number"],
            main_product_index=product_data["main_product_index"],
            product_molecular_weights=product_data["product_molecular_weights"],
            product_masses=product_data["product_masses"],
            rounded_product_masses=product_data["rounded_product_masses"],
            ae=sustainability_data["ae"],
            ae_flag=sustainability_data["ae_flag"],
            element_sustainability=sustainability_data["element_sustainability"],
            element_sustainability_flag=sustainability_data[
                "element_sustainability_flag"
            ],
            reactant_hazard_sentences=reactant_data["reactant_hazard_sentences"],
            reactant_hazard_ratings=reactant_data["reactant_hazard_ratings"],
            reactant_hazard_colors=reactant_data["reactant_hazard_colours"],
            reactant_risk_colors=reactant_data["reactant_risk_colours"],
            reactant_exposure_potentials=reactant_data["reactant_exposure_potentials"],
            reactant_risk_ratings=reactant_data["reactant_risk_ratings"],
            reagent_hazard_sentences=reagent_data["reagent_hazard_sentences"],
            reagent_hazard_ratings=reagent_data["reagent_hazard_ratings"],
            reagent_hazard_colors=reagent_data["reagent_hazard_colours"],
            reagent_risk_colors=reagent_data["reagent_risk_colours"],
            reagent_exposure_potentials=reagent_data["reagent_exposure_potentials"],
            reagent_risk_ratings=reagent_data["reagent_risk_ratings"],
            solvent_primary_keys=solvent_data["solvent_primary_keys_str"],
            solvent_hazard_sentences=solvent_data["solvent_hazard_sentences"],
            solvent_hazard_ratings=solvent_data["solvent_hazard_ratings"],
            solvent_exposure_potentials=solvent_data["solvent_exposure_potentials"],
            solvent_risk_ratings=solvent_data["solvent_risk_ratings"],
            solvent_hazard_colors=solvent_data["solvent_hazard_colours"],
            solvent_risk_colors=solvent_data["solvent_risk_colours"],
            product_hazard_sentences=product_data["product_hazard_sentences"],
            product_hazard_ratings=product_data["product_hazard_ratings"],
            product_exposure_potentials=product_data["product_exposure_potentials"],
            product_risk_ratings=product_data["product_risk_ratings"],
            product_hazard_colors=product_data["product_hazard_colours"],
            product_risk_colors=product_data["product_risk_colours"],
            risk_rating=risk_data["risk_rating"],
            risk_color=risk_data["risk_colour"],
            number_of_solvents=solvent_data["number_of_solvents"],
        )
        return jsonify({"summary": summary_table})
    else:
        pass
    return jsonify(
        {"summary": "Ensure you have entered all the necessary information!"}
    )  # otherwise it shows this message


@summary_bp.route("/element_sustainability", methods=["POST", "GET"])
@login_required
def element_sustainability() -> Response:
    # must be logged in
    return render_template("element_sustainability.html")


@summary_bp.route("/pdf", methods=["POST", "GET"])
@login_required
def pdf():
    """
    Saves the autogenerated PDF for a reaction, overwriting the old one if present.
    """
    services.auth.reaction_files(permission_level="edit")
    new_upload = services.file_attachments.UploadExperimentDataFiles(
        request, autogenerated_file=True
    )
    new_upload.validate_files()
    new_upload.remove_duplicate_autogenerated_summaries()
    new_upload.remove_duplicate_autogenerated_summaries()
    new_upload.save_validated_files()
    return "", 204


@summary_bp.route("/get_file_attachment_list", methods=["POST"])
@login_required
def get_file_attachment_list():
    """
    Gets a list of file attachments for the current reaction identified from the request data
    """
    services.auth.reaction_files(permission_level="view")
    reaction = services.reaction.get_current_from_request()
    file_attachments = sorted(
        reaction.file_attachments, key=lambda x: not x.autogenerated
    )
    file_attachments_dict_list = [
        {"name": x.display_name, "uuid": x.uuid, "autogenerated": x.autogenerated}
        for x in file_attachments
    ]
    return jsonify({"file_attachments": file_attachments_dict_list})
