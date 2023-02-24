from flask_login import login_required
from sources.export_data import export_data_bp
from flask import render_template, flash, redirect, url_for
from pony.orm import select, desc
from sources import db
from sources.auxiliary import get_workgroups, get_notification_number, security_member_workgroup_workbook
from flask import send_file
import pandas as pd
import csv
import io
import ast


@export_data_bp.route('/export_data_csv/<workgroup>/<workbook>', methods=['GET', 'POST'])
@login_required
def export_data_csv(workgroup, workbook):
    # get all reaction data for user
    reaction_list = select(reaction for reaction in db.Reaction if reaction.workbooks.name == workbook and
                           reaction.workbooks.group.name == workgroup and reaction.status == "active")[:]
    # convert to list of dictionaries
    reaction_list = [x.to_dict() for x in reaction_list]
    # remove empty columns
    columns = ['workbooks', 'reaction_class', 'date_reaction', 'green_metric', 'status']
    for reaction in reaction_list:
        for col in columns:
            del reaction[col]
    # convert reaction table and summary table to dictionaries and remove summary to print
    for reaction in reaction_list:
        reaction["reaction_table_data"] = ast.literal_eval(reaction["reaction_table_data"])
        reaction["summary_table_data"] = ast.literal_eval(reaction["summary_table_data"])
        try:
            del reaction["summary_table_data"]["summary_to_print"]
        except KeyError:
            pass
    # make data frame
    df = pd.DataFrame(reaction_list)
    # get creator usernames
    creators = []
    for creator in df["creator"]:
        if db.Person[creator].user:
            creators.append(db.Person[creator].user.username)
    df["creator"] = creators
    # reactants, reagents, products, or solvents change to empty

    def clear_empty(column):
        ls_column = df[column].to_list()
        ls_column = [str(x).replace(", ''", "") for x in ls_column]
        ls_column = [str(x).replace("''", "") for x in ls_column]
        ls_column = [str(x).replace("[]", "") for x in ls_column]
        return ls_column

    df["reactants"] = clear_empty("reactants")
    df["products"] = clear_empty("products")
    df["reagents"] = clear_empty("reagents")
    df["solvent"] = clear_empty("solvent")
    # make selects empty
    df = df.replace(to_replace="\"-select-\"", value="\"\"", regex=True)

    # get solvent names
    reaction_table_list = df["reaction_table_data"].tolist()
    solvent_list = [x["solvent_names"] for x in reaction_table_list]
    df["solvent"] = solvent_list

    # get solvent sustainability
    flags = {4: "Recommended", 3: "Problematic", 2: "Hazardous", 1: "Highly Hazardous", 5: 'non-chem21'}
    solvent_hazards = []
    for sol_list in df["solvent"].tolist():
        if sol_list != ['']:
            solvent_hazards.append([flags[select(x.flag for x in db.Solvent if x.name == n).first()] for n in sol_list])
        else:
            solvent_hazards.append("")
    df.insert(11, "solvent sustainability", solvent_hazards)
    # remove none values
    df = df.replace("['None']", "")
    # replace indices in reaction and summary table for actual values
    replace_elements = {"3": "+500 years", "2": "50-500 years", "1": "5-50 years", "0": ""}
    try:
        for summary in df["summary_table_data"].tolist():
            summary["element_sustainability"] = replace_elements[summary["element_sustainability"]]

        replace_isolation = {"0": "",
                             "1": "Column",
                             "2": "HPLC",
                             "3": "Ion exchange",
                             "4": "Crystallization",
                             "5": "Filtration",
                             "6": "Multiple recryst.",
                             "7": "Distillation < 140 degC",
                             "8": "Distillation > 140 degC"}
        for summary in df["summary_table_data"].tolist():
            summary["isolation_method"] = replace_isolation[summary["isolation_method"]]

        # replace all physical form indices
        replace_phys_forms = {"0": "", "1": "Dense Solid", "2": "Non-volatile liquid (b.p. > 130 degC)", "3": "Unknown",
                              "4": "Dusty Solid", "5": "Lyophilised solid", "6": "Volatile liquid (b.p. 70-130 degC)",
                              "7": "Gas", "8": "Highly volatile liquid (b.p. < 70 degC)", "9": "Aerosol",
                              "10": "Solution that promotes skin absorption"}
        phys_form_substances = ["reactant_physical_forms", "reagent_physical_forms", "solvent_physical_forms",
                                "product_physical_forms"]
        for reaction in df["reaction_table_data"].tolist():
            for subs in phys_form_substances:
                reaction[subs] = [replace_phys_forms[x] for x in reaction[subs]]

        # reformat extra reaction and summary data
        to_replace_replace = [["\",\"", "\":\""],
                              ["\"key\":", ""],
                              ["\"value\":", ""],
                              ["{", ""],
                              ["}", ""],
                              ["[", ""],
                              ["]", ""],
                              ["hazard\":\"1\"", "hazard\":\"1-Slight\""],
                              ["hazard\":\"2\"", "hazard\":\"2-Serious\""],
                              ["hazard\":\"3\"", "hazard\":\"3-Major\""],
                              ["risk\":\"1\"", "risk\":\"1-Low likelihood\""],
                              ["risk\":\"2\"", "risk\":\"2-Possible\""],
                              ["risk\":\"3\"", "risk\":\"3-Frequent Occur\""],
                              ["consequences\":\"1\"", "consequences\":\"1-Individual\""],
                              ["consequences\":\"2\"", "consequences\":\"2-Local Labs\""],
                              ["consequences\":\"3\"", "consequences\":\"3-Building Wide\""],
                              ["hazard-acceptable", "Recommended"],
                              ["hazard-warning", "Problematic"],
                              ["hazard-hazardous", "Hazardous"],
                              ["hazard-reset-hazard", ""],
                              [" to-export", ""]]
        df["summary_table_data_supplementary"] = [x["to_export"] for x in df["summary_table_data"].tolist()]
        print(df["summary_table_data_supplementary"].tolist())
        for t_r in to_replace_replace:
            df["summary_table_data_supplementary"] = [x.replace(t_r[0], t_r[1]) for x in
                                                      df["summary_table_data_supplementary"].tolist()]
        df["summary_table_data_supplementary"] = ["{" + x for x in df["summary_table_data_supplementary"].tolist()]
        df["summary_table_data_supplementary"] = [x + "}" for x in df["summary_table_data_supplementary"].tolist()]
        # remove to_export raw data
        for reaction in reaction_list:
            del reaction["summary_table_data"]["to_export"]
        # put summary table data together in the same column
        df["summary_table_data"] = [str(df["summary_table_data"].tolist()[i]) +
                                    df["summary_table_data_supplementary"].tolist()[i] for i in range(df.shape[0])]
        df["summary_table_data"] = [x.replace("}{", ",") for x in df["summary_table_data"].tolist()]
        df["summary_table_data"] = [x.replace("\"", "'") for x in df["summary_table_data"].tolist()]
        df = df.drop(columns=["summary_table_data_supplementary"])
    except KeyError:
        pass

    column_names = list(df)
    df = df.values.tolist()
    proxy = io.StringIO()

    writer = csv.writer(proxy)
    writer.writerow(column_names)
    for row in df:
        writer.writerow(row)

    # Creating the byteIO object from the StringIO Object
    mem = io.BytesIO()
    mem.write(proxy.getvalue().encode())
    # seeking was necessary
    mem.seek(0)
    proxy.close()

    return send_file(
        mem,
        as_attachment=True,
        download_name='reaction_data.csv',
        mimetype='text/csv'
    )


@export_data_bp.route('/export_data_pdf/<workgroup>/<workbook>/<sort_crit>', methods=['GET', 'POST'])
@export_data_bp.route('/export_data_pdf/<workgroup>/<workbook>/<sort_crit>/<reaction_name>', methods=['GET', 'POST'])
@login_required
def export_data_pdf(workgroup, workbook, sort_crit, reaction_name=None):
    # must be logged in and member of workbook and workgroup
    if not security_member_workgroup_workbook(workgroup, workbook):
        flash("You do not have permission to view this page")
        return redirect(url_for('main.index'))
    if sort_crit == "time":
        reaction_list = select(reaction for reaction in db.Reaction if reaction.workbooks.name == workbook and
                               reaction.workbooks.group.name == workgroup and reaction.status == "active").order_by(desc(db.Reaction.time_of_creation))[:]
    elif sort_crit == "AZ":
        reaction_list = select(reaction for reaction in db.Reaction if reaction.workbooks.name == workbook and
                               reaction.workbooks.group.name == workgroup and reaction.status == "active").order_by(db.Reaction.name)[:]
    else:
        reaction_list = select(reaction for reaction in db.Reaction if reaction.workbooks.name == workbook and
                               reaction.workbooks.group.name == workgroup and reaction.name == reaction_name)[:]
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    return render_template("_summary_print.html", reaction_list=reaction_list, sort_crit=sort_crit,
                           workgroups=workgroups, notification_number=notification_number)

