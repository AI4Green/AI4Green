#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module receives data from the reaction table
via POST request and renders the summary template
"""
from flask import render_template, request, jsonify
import json
# render_template renders html templates
# request parses incoming request data and gives access to it
# jsonify is used to send a JSON response to the browser
from sources.summary import summary_bp  # imports the blueprint of the summary route
from sources import auxiliary, db  # imports the module with auxiliary functions
from pony.orm import select
import ast
import rdkit
from flask_login import login_required


# Processing data from the reaction table and creating summary with H&S report
@summary_bp.route('/_summary', methods=['POST'])
@login_required
def summary():
    # must be logged in
    """This function receives the reaction information from browser, calculates
    green metrics, gives hazard information, and renders the summary"""
    summary_table_data = str(request.form['js_summary_table_data'])
    if summary_table_data != "no data":
        summary_table_data = ast.literal_eval(summary_table_data)
    mass_factor = {'g': 1, 'mg': 10 ** (-3), 'μg': 10 ** (-6)}  # mass factor dictionary to convert mass units
    category_rate = {'': 0, 'L': 1, 'M': 2, 'H': 3, 'VH': 4}  # hazard category rate dictionary
    # Gets selected units from the reaction table
    amount_unit = str(request.form['amountUnit'])  # mol, mmol   -
    volume_unit = str(request.form['volumeUnit'])  # L, mL etc    | - These three are the units used for ALL reactants + reagents
    mass_unit = str(request.form['massUnit'])  # g, mg, etc       -
    solvent_volume_unit = str(request.form['solventVolumeUnit'])  # L, ml, etc
    product_mass_unit = str(request.form['productMassUnit'])  # g, mg, etc
    reactant_mass_factor = mass_factor[mass_unit]
    product_mass_factor = mass_factor[product_mass_unit]

    # Gets reactant data from the reaction table
    reactants = auxiliary.get_data('reactants')  # reactants='Ethanol;Methanol;...'
    reactant_molecular_weights = auxiliary.get_data('reactantMolecularWeights')  # reactantMolecularWeights='4.56;7.86'
    reactant_densities = auxiliary.get_data('reactantDensities')  # reactantDensities='9.8;890.02'
    reactant_concentrations = auxiliary.get_data('reactantConcentrations')
    reactant_equivalents = auxiliary.get_data('reactantEquivalents')
    reactant_amounts = auxiliary.get_data('reactantAmounts')
    rounded_reactant_amounts = auxiliary.get_data('roundedReactantAmounts')
    reactant_volumes = auxiliary.get_data('reactantVolumes')
    rounded_reactant_volumes = auxiliary.get_data('roundedReactantVolumes')
    reactant_masses = auxiliary.get_data('reactantMasses')
    rounded_reactant_masses = auxiliary.get_data('roundedReactantMasses')
    reactant_hazards = auxiliary.get_data('reactantHazards')  # reactantHazards='H344;H342'
    reactant_physical_forms = auxiliary.get_data('reactantPhysicalForms')
    reactant_mass_sum = float(request.form['reactantMassSum'])
    reactant_molecular_weight_sum = float(request.form['reactantMolecularWeightSum'])
    reactant_primary_keys = auxiliary.get_data('reactantPrimaryKeys')
    reactant_primary_keys = ", ".join(reactant_primary_keys)

    # Gets reagent data from the reaction table
    reagent_table_numbers = auxiliary.get_data('reagentTableNumbers')
    reagents = auxiliary.get_data('reagents')
    reagent_molecular_weights = auxiliary.get_data('reagentMolecularWeights')
    reagent_densities = auxiliary.get_data('reagentDensities')
    reagent_concentrations = auxiliary.get_data('reagentConcentrations')
    reagent_equivalents = auxiliary.get_data('reagentEquivalents')
    reagent_amounts = auxiliary.get_data('reagentAmounts')
    rounded_reagent_amounts = auxiliary.get_data('roundedReagentAmounts')
    reagent_volumes = auxiliary.get_data('reagentVolumes')
    rounded_reagent_volumes = auxiliary.get_data('roundedReagentVolumes')
    reagent_masses = auxiliary.get_data('reagentMasses')
    rounded_reagent_masses = auxiliary.get_data('roundedReagentMasses')
    reagent_hazards = auxiliary.get_data('reagentHazards')
    reagent_physical_forms = auxiliary.get_data('reagentPhysicalForms')
    reagent_mass_sum = float(request.form['reagentMassSum'])
    reagent_molecular_weight_sum = float(request.form['reagentMolecularWeightSum'])
    reagent_primary_keys_ls = auxiliary.get_data('reagentPrimaryKeys')
    if reagent_primary_keys_ls == ['']:
        reagent_primary_keys_ls = ['0']
    reagent_primary_keys = ", ".join(reagent_primary_keys_ls)

    # Gets solvent data from the reaction table
    solvent_table_numbers = auxiliary.get_data('solventTableNumbers')
    solvents = auxiliary.get_data('solvents')
    solvent_volumes = auxiliary.get_data('solventVolumes')
    solvent_hazards = auxiliary.get_data('solventHazards')
    solvent_physical_forms = auxiliary.get_data('solventPhysicalForms')
    number_of_solvents = request.form['numberOfSolvents']
    solvent_primary_keys_ls = auxiliary.get_data('solventPrimaryKeys')
    if solvent_primary_keys_ls == ['']:
        solvent_primary_keys_ls = ['0']
    solvent_primary_keys = ", ".join(solvent_primary_keys_ls)

    # Gets product data from the reaction table
    product_table_numbers = list(filter(None, auxiliary.get_data('productTableNumbers')))
    product_table_numbers = [int(x) for x in product_table_numbers]
    products = auxiliary.get_data('products')
    product_masses = auxiliary.get_data('productMasses')
    rounded_product_masses = auxiliary.get_data('roundedProductMasses')
    product_molecular_weights = auxiliary.get_data('productMolecularWeights')
    product_hazards = auxiliary.get_data('productHazards')
    product_physical_forms = auxiliary.get_data('productPhysicalForms')
    product_primary_keys = auxiliary.get_data('productPrimaryKeys')
    product_primary_keys = ", ".join(product_primary_keys)

    # Gets main product from the reaction table
    main_product_table_number = int(request.form['mainProductTableNumber'])
    main_prod_idx = 'unassigned'
    for idx, product_table_num in enumerate(product_table_numbers):
        if product_table_num == main_product_table_number:
            main_prod_idx = int(idx)
            break
    # check all the requirement information has been typed into the reaction table
    # this is limiting reagent mass, any other reagent equivalents, physical forms and solvent volume
    for mass in reactant_masses:
        if mass == "" or 0:
            print("1")
            return jsonify({'summary': 'Ensure you have entered all the necessary information!'})
    for equi in reactant_equivalents:
        if equi == "":
            print("2")
            return jsonify({'summary': 'Ensure you have entered all the necessary information!'})
    if reagent_table_numbers[0]:
        for equi in reagent_equivalents:
            if equi == "":
                print("3")
                return jsonify({'summary': 'Ensure you have entered all the necessary information!'})
    if reagent_table_numbers[0]:
        for name in reagents:
            if name == "":
                print("4")
                return jsonify({'summary': 'Ensure you have entered all the necessary information!'})
        # if name not in database, no hazard code
        for haz in reagent_hazards:
            if haz == "":
                print("5")
                return jsonify({'summary': 'Ensure you have entered all the necessary information!'})

    for phys_form in reactant_physical_forms:
        if phys_form == "-select-":
            print("6")
            return jsonify({'summary': 'Ensure you have entered all the necessary information!'})
    for phys_form in reagent_physical_forms:
        if phys_form == "-select-":
            print("7")
            return jsonify({'summary': 'Ensure you have entered all the necessary information!'})
    for phys_form in solvent_physical_forms:
        if phys_form == "-select-":
            print("8")
            return jsonify({'summary': 'Ensure you have entered all the necessary information!'})
    for phys_form in product_physical_forms:
        if phys_form == "-select-":
            print("9")
            return jsonify({'summary': 'Ensure you have entered all the necessary information!'})
    if solvent_table_numbers[0]:
        for vol in solvent_volumes:
            if vol == "" or 0:
                print("10")
                return jsonify({'summary': 'Ensure you have entered all the necessary information!'})
    if solvent_table_numbers[0]:
        for sol in solvents:
            if sol == "" or 0:
                print("11")
                return jsonify({'summary': 'Ensure you have entered all the necessary information!'})

    # Calculates green metrics
    pmi = round((reactant_mass_sum + reagent_mass_sum) * reactant_mass_factor /
                (float(product_masses[main_prod_idx]) * product_mass_factor), 1) \
        if float(product_masses[main_prod_idx]) > 0 else 0  # product mass intensity
    ae = round(100 * float(product_molecular_weights[main_prod_idx]) /
               (reactant_molecular_weight_sum + reagent_molecular_weight_sum), 1) \
        if (reactant_molecular_weight_sum + reagent_molecular_weight_sum) > 0 else 0  # atom economy
    ae_flag = auxiliary.metric_flag(ae)  # atom economy color flag
    # rme = round(100 * (float(product_masses[main_prod_idx]) * product_mass_factor) / ((reactant_mass_sum + reagent_mass_sum) * reactant_mass_factor), 1) \
    #    if (reactant_mass_sum + reagent_mass_sum) > 0 else 0  # reaction mass efficiency
    # rme_flag = auxiliary.metric_flag(rme)  # reaction mass efficiency color flag
    # oe = round(100 * rme / ae, 1) if ae > 0 else 0  # optimum efficiency
    # ef = round(pmi - 1, 1) if pmi > 1 else 0  # environmental factor (E-factor)
    # mp = round(100 / pmi, 1) if pmi > 0 else 0  # mass productivity
    # get reaction smiles and convert to a list of smiles involved in reaction
    reaction_smiles = str(request.form['reactionSmiles']).split(' |')[0]
    reaction_smiles_ls = reaction_smiles.replace('>>', '.').split('.')
    # get reagent and solvent smiles
    reagent_smiles_ls = [select(y.smiles for y in db.Compound if y.id == x).first() for x in reagent_primary_keys_ls]
    solvent_smiles_ls = [select(y.smiles for y in db.Compound if y.id == x).first() for x in solvent_primary_keys_ls]
    full_reaction_smiles_ls = reaction_smiles_ls + reagent_smiles_ls + solvent_smiles_ls
    full_reaction_smiles_ls1 = [x for x in full_reaction_smiles_ls if x is not None]
    # get list of elements in reaction from the smiles list by converting to atoms and then chemical symbols
    element_symbols = set()
    for component in full_reaction_smiles_ls1:
        mol = rdkit.Chem.MolFromSmiles(component)
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            element_symbols.add(symbol)
    # get sustainability colour and set element_sustainability variable according to least sustainable element present
    element_sustainability_set = set([select(y.colour for y in db.Element if y.symbol == symbol).first()
                                      for symbol in element_symbols])
    element_sustainability = '-select-'
    element_sustainability_flag = 'hazard-reset-hazard'
    if 'red' in element_sustainability_set:
        element_sustainability = '5-50 years'
        element_sustainability_flag = 'hazard-hazardous'
    elif 'yellow' in element_sustainability_set:
        element_sustainability = '50-500 years'
        element_sustainability_flag = 'hazard-warning'
    elif 'lime' in element_sustainability_set:
        element_sustainability = '+500 years'
        element_sustainability_flag = 'hazard-acceptable'
    # Hazard summary
    # Reactant hazards
    reactant_hazard_sentences, reactant_hazard_ratings, \
    reactant_hazard_colors, reactant_risk_colors, \
    reactant_exposure_potentials, reactant_risk_ratings, reactant_total_rate \
        = auxiliary.compound_hazard_data(reactants, reactant_hazards, reactant_physical_forms)

    # Reagent hazards
    reagent_hazard_sentences, reagent_hazard_ratings, \
    reagent_hazard_colors, reagent_risk_colors, \
    reagent_exposure_potentials, reagent_risk_ratings, reagent_total_rate \
        = auxiliary.compound_hazard_data(reagents, reagent_hazards, reagent_physical_forms)
    # solvent hazards
    solvent_hazard_sentences, solvent_hazard_ratings, \
    solvent_hazard_colors, solvent_risk_colors, \
    solvent_exposure_potentials, solvent_risk_ratings, solvent_total_rate \
        = auxiliary.compound_hazard_data(solvents, solvent_hazards, solvent_physical_forms)

    # Product hazard
    product_hazard_sentences, product_hazard_ratings, \
    product_hazard_colors, product_risk_colors, \
    product_exposure_potentials, product_risk_ratings, product_total_rate \
        = auxiliary.compound_hazard_data(products, product_hazards, product_physical_forms)
    """"""

    total_rate = reactant_total_rate + reagent_total_rate + solvent_total_rate + product_total_rate
    max_total_rate = int(max(total_rate))  # max total hazard rate
    risk_rating = list(category_rate.keys())[max_total_rate]  # resulting hazard rating
    risk_color = 'hazard-hazardous' if risk_rating == 'VH' else 'hazard-reset-hazard'  # color code for the hazard rating

    # Solvent flags
    flag_rate = {1: 'hazard-highly-hazardous', 4: 'hazard-acceptable', 3: 'hazard-warning', 2: 'hazard-hazardous',
                 5: 'non-chem21'}  # flag rate dictionary
    solvent_flags = []  # solvent flag list
    if solvents[0]:
        for solvent in solvents:
            solvent_flag = select(x.flag for x in db.Solvent if x.name == solvent).first()
            solvent_flags.append(flag_rate[solvent_flag])  # appends solvent flag to their list

    if number_of_solvents == "0":  # if no solvents have been chosen
        number_of_solvents = 1  # then it shows only one empty cell
        solvents = [' ']

    # if product mass and reactant mass sum are calculated, then it forms a summary table
    if product_masses and reactant_mass_sum:
        summary_table = render_template("_summary_table.html",
                                        amount_unit=amount_unit, volume_unit=volume_unit,
                                        mass_unit=mass_unit, solvent_volume_unit=solvent_volume_unit,
                                        product_mass_unit=product_mass_unit, reactants=reactants,
                                        reactant_primary_keys=reactant_primary_keys,
                                        reagent_primary_keys=reagent_primary_keys,
                                        reagents=reagents, reagent_table_numbers=reagent_table_numbers,
                                        reagent_molecular_weights=reagent_molecular_weights,
                                        reagent_densities=reagent_densities,
                                        reagent_concentrations=reagent_concentrations,
                                        reagent_equivalents=reagent_equivalents, reagent_hazards=reagent_hazards,
                                        reagent_amounts=reagent_amounts,
                                        rounded_reagent_amounts=rounded_reagent_amounts,
                                        reagent_volumes=reagent_volumes,
                                        rounded_reagent_volumes=rounded_reagent_volumes,
                                        reagent_masses=reagent_masses, rounded_reagent_masses=rounded_reagent_masses,
                                        solvents=solvents, solvent_volumes=solvent_volumes,
                                        solvent_table_numbers=solvent_table_numbers,
                                        solvent_flags=solvent_flags,
                                        products=products, product_table_numbers=product_table_numbers,
                                        reactant_molecular_weights=reactant_molecular_weights,
                                        reactant_densities=reactant_densities,
                                        reactant_concentrations=reactant_concentrations,
                                        reactant_equivalents=reactant_equivalents,
                                        reactant_amounts=reactant_amounts,
                                        rounded_reactant_amounts=rounded_reactant_amounts,
                                        reactant_volumes=reactant_volumes,
                                        rounded_reactant_volumes=rounded_reactant_volumes,
                                        reactant_masses=reactant_masses,
                                        rounded_reactant_masses=rounded_reactant_masses,
                                        product_primary_keys=product_primary_keys,
                                        main_product_table_number=main_product_table_number,
                                        main_product_index=main_prod_idx,
                                        product_molecular_weights=product_molecular_weights,
                                        product_masses=product_masses,
                                        rounded_product_masses=rounded_product_masses,
                                        ae=ae, ae_flag=ae_flag,  # rme=rme, rme_flag=rme_flag,
                                        # oe=oe, pmi=pmi, ef=ef, mp=mp,
                                        element_sustainability=element_sustainability,
                                        element_sustainability_flag=element_sustainability_flag,
                                        reactant_hazard_sentences=reactant_hazard_sentences,
                                        reactant_hazard_ratings=reactant_hazard_ratings,
                                        reactant_hazard_colors=reactant_hazard_colors,
                                        reactant_risk_colors=reactant_risk_colors,
                                        reactant_exposure_potentials=reactant_exposure_potentials,
                                        reactant_risk_ratings=reactant_risk_ratings,
                                        reagent_hazard_sentences=reagent_hazard_sentences,
                                        reagent_hazard_ratings=reagent_hazard_ratings,
                                        reagent_hazard_colors=reagent_hazard_colors,
                                        reagent_risk_colors=reagent_risk_colors,
                                        reagent_exposure_potentials=reagent_exposure_potentials,
                                        reagent_risk_ratings=reagent_risk_ratings,
                                        solvent_primary_keys=solvent_primary_keys,
                                        solvent_hazard_sentences=solvent_hazard_sentences,
                                        solvent_hazard_ratings=solvent_hazard_ratings,
                                        solvent_exposure_potentials=solvent_exposure_potentials,
                                        solvent_risk_ratings=solvent_risk_ratings,
                                        solvent_hazard_colors=solvent_hazard_colors,
                                        solvent_risk_colors=solvent_risk_colors,
                                        product_hazard_sentences=product_hazard_sentences,
                                        product_hazard_ratings=product_hazard_ratings,
                                        product_exposure_potentials=product_exposure_potentials,
                                        product_risk_ratings=product_risk_ratings,
                                        product_hazard_colors=product_hazard_colors,
                                        product_risk_colors=product_risk_colors,
                                        risk_rating=risk_rating, risk_color=risk_color,
                                        number_of_solvents=number_of_solvents,
                                        summary_table_data=json.dumps(summary_table_data))
        return jsonify({'summary': summary_table})
    else:
        pass
    return jsonify({
        'summary': 'Ensure you have entered all the necessary information!'})  # otherwise it shows this message


@summary_bp.route('/element_sustainability', methods=['POST', 'GET'])
@login_required
def element_sustainability():
    # must be logged in
    return render_template('element_sustainability.html')
