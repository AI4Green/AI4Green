"""
The code here is for constructing the chemical database using the XML from PubChem
"""

import contextlib
import csv
import re
import sqlite3 as sql
from datetime import datetime
from statistics import median, mode
from typing import Union

import pytz
from chemspipy import ChemSpider
from compound_database.auxiliary import compound_database_dir
from lxml import etree as et
import os
from rdkit import Chem
from sources import models
from sources.extensions import db
from utilities import read_yaml, basedir

ns = {
    "xmlns": "http://pubchem.ncbi.nlm.nih.gov/pug_view",  # Defines the url used in the xml tags
    "xs": "http://pubchem.ncbi.nlm.nih.gov/pug_view",
}
current_time = datetime.now(pytz.timezone("Europe/London")).replace(tzinfo=None)


def get_cas_num(root):
    """Returns the CAS number when given the root for a single record.
    :param root:
    :return cas_number, cas_number_ls:
    """

    for sections in root.findall("xmlns:Section/xmlns:Section", ns):
        # Finds where the CAS number is stored in the XML
        for headings in sections.findall("xmlns:TOCHeading", ns):
            if headings.text == "CAS":
                # Saves the CAS numbers as a list
                cas_number_ls = [
                    casNum.text
                    for casNum in sections.findall(
                        "xmlns:Information/xmlns"
                        ":Value/xmlns:StringWithMarkup/xmlns:String",
                        ns,
                    )
                ]
                # Finds the mode of the cas numbers.
                cas_number = mode(cas_number_ls)
                return cas_number
            else:
                return None  # returning both as none means the record will be skipped if no cas is found


def cas_verification(
    cas_number: str,
    model: Union[models.Compound, models.UpdateCompound] = models.Compound,
) -> bool:
    """
    Determine if a CAS number is valid or not.

    Remarks:
    A CAS number is invalid if it is not present, incorrect length, or already present in the database.

    Args:
        cas_number: CAS number of a compound
        model: The model to check against - this reflects the database bind used.

    returns:
        bool: Whether the CAS number is valid.

    """
    duplicate_cas_check = (
        db.session.query(model).filter(model.cas == cas_number).first()
    )
    if duplicate_cas_check is not None:
        return False

    if cas_number is None:
        return False
    elif not 7 <= len(cas_number) <= 13:
        print("wrong length cas - compound omitted from database")
        return False
    else:
        return True


def get_compound_name(root):
    """This function retrieves the 'main' name for a compound that is in the record title of the xml file
    :param root:
    :return compound_name:

    """

    for tags in root.findall("xmlns:RecordTitle", ns):
        return tags.text


def get_cid(root):
    """This function retrieves the pub chem compound identifier (CID) which is unique to a compound.
    Type should be either RecordNumber or RecordTitle
    :param root:
    :return record_info:
    """
    for tags in root.findall("xmlns:RecordNumber", ns):
        return tags.text


def get_identifier(root, id_type):
    """This function is used to retrieve the InChI and InChI key from the xml file
    :param root:
    :param id_type:
    :return id:
    """
    for sections in root.findall("xmlns:Section/xmlns:Section", ns):
        for headings in sections.findall("xmlns:TOCHeading", ns):
            if headings.text == id_type:
                for id_values in sections.findall(
                    "xmlns:Information/xmlns:Value/xmlns:StringWithMarkup/xmlns:String",
                    ns,
                ):
                    return id_values.text


def get_smiles(inchi_str):
    """
    Function for converting the InChI string into a Smiles string
    :param inchi_str:
    :return:
    """

    with contextlib.suppress(Exception):
        mol = Chem.inchi.MolFromInchi(
            inchi_str,
            sanitize=True,
            removeHs=True,
            logLevel=None,
            treatWarningAsError=False,
        )
        return Chem.MolToSmiles(mol)
    with contextlib.suppress(Exception):
        mol = Chem.inchi.MolFromInchi(
            inchi_str, sanitize=False, logLevel=None, treatWarningAsError=False
        )
        return Chem.MolToSmiles(mol)
    try:
        cs = ChemSpider("584QozVSGkkJn8YulyGGSthFGIgEB6PJ")
        return cs.convert(input=inchi_str, input_format="InChI", output_format="SMILES")
    except Exception:
        return "smiles error"


def get_mol_weight(root):
    """This function returns the molecular weight when given the root for a single record
    :param root:
    :return molecular_weight:
    """

    for sections in root.findall("xmlns:Section", ns):
        # Finds where the molecular weight is stored in the xml
        for headings in sections.findall("xmlns:TOCHeading", ns):
            if headings.text == "Molecular Weight":
                # Saves the molecular weight and converts it to a string.
                molecular_weight_ls = [
                    mol_weight.text
                    for mol_weight in sections.findall(
                        "xmlns:Information/xmlns:Value/xmlns:StringWithMarkup/"
                        "xmlns:String",
                        ns,
                    )
                ]
                molecular_weight = molecular_weight_ls[0]
                return float(molecular_weight)


def find_echa(root):
    """This function finds the ECHA reference number used for the hazard data.
    This means all hazard data extracted is the data supplied to PubChem by the ECHA
    :param root:
    :return ref_num:
    """

    for sections in root.findall("xmlns:Section/xmlns:Information", ns):
        for references in sections.findall("xmlns:Name", ns):
            if references.text == "ECHA C&L Notifications Summary":
                ref_num_ls = [
                    ref_nums.text
                    for ref_nums in sections.findall("xmlns:ReferenceNumber", ns)
                ]
                return ref_num_ls[0]
        # alternative name
        for references in sections.findall("xmlns:Name", ns):
            if references.text == "Signal":
                ref_num_ls = [
                    ref_nums.text
                    for ref_nums in sections.findall("xmlns:ReferenceNumber", ns)
                ]
                return ref_num_ls[0]


def get_hazards(root, echa_ref):
    """Returns hazard codes when given the root for a single record and the reference number for the echa data
    :param root:
    :param echa_ref:
    :return hazard_codes:
    """

    for sections in root.findall("xmlns:Section/xmlns:Information", ns):
        # Finds where the hazard codes are stored
        for headings in sections.findall("xmlns:Name", ns):
            if headings.text == "GHS Hazard Statements":
                # Finds the ECHA supplied data using the reference number and saves them as a list
                for hazard_suppliers in sections.findall("xmlns:ReferenceNumber", ns):
                    if hazard_suppliers.text == echa_ref:
                        hazard_codes_ls = [
                            haz_codes.text
                            for haz_codes in sections.findall(
                                "xmlns:Value/xmlns:StringWithMarkup/xmlns:String", ns
                            )
                        ]

                        # Formats the hazard codes by separating the different codes then only saving the H codes.
                        hazard_codes = ""

                        for hazards in hazard_codes_ls:
                            hazards_text = hazards.split(" ")
                            ind_hazard_codes = hazards_text[0]
                            ind_hazard_codes = ind_hazard_codes.split(":")[0]
                            ind_hazard_codes += "-"
                            hazard_codes += ind_hazard_codes
                        return hazard_codes[:-1]


def get_mol_formula(root):
    """
    Extracts the molecular formula for a single record
    :param root:
    :return mol_formula:
    """

    for sections in root.findall("xmlns:Section", ns):
        for headings in sections.findall("xmlns:TOCHeading", ns):
            if headings.text == "Molecular Formula":
                for values in sections.findall(
                    "xmlns:Information/xmlns:Value/xmlns:StringWithMarkup/xmlns:String",
                    ns,
                ):
                    return values.text


def get_phys_prop(root, property_input):
    """This function is used to retrieve the physical properties from the xml file (temperatures, form and density)
    :param root:
    :param property_input:
    :return property:
    """

    for sections in root.findall("xmlns:Section/xmlns:Section", ns):
        for headings in sections.findall("xmlns:TOCHeading", ns):
            if headings.text == property_input:
                for entries in sections.findall(
                    "xmlns:Information/xmlns:Description", ns
                ):
                    if entries.text == "PEER REVIEWED":
                        phys_properties_ls = [
                            values.text
                            for values in sections.findall(
                                "xmlns:Information/xmlns"
                                ":Value/xmlns"
                                ":StringWithMarkup/xmlns"
                                ":String",
                                ns,
                            )
                        ]

                        temperatures = {
                            "Boiling Point",
                            "Melting Point",
                            "Flash Point",
                            "Autoignition Temperature",
                        }

                        if phys_properties_ls:
                            if property_input in temperatures:
                                phys_property = temp_eval(phys_properties_ls)
                                return phys_property

                            elif property_input == "Density":
                                phys_property = density_eval(phys_properties_ls)
                                return phys_property


def temp_eval(lls):
    """Converts any temperatures in fahrenheit and checks the temperatures are consistent with one another
    :param lls:
    :return rounded temp:
    """

    def extract_temp(lls):
        """Finds the temperature value, identifies when temperature is presented as a range and
        separates temperatures into celsius and fahrenheit.

        :param lls:
        :return celsius_temps, fahrenheit_temps:
        """
        c_temp = []
        f_temp = []
        decomposition = ["decomposes", "Decomposes", "decomposition", "Decomposition"]
        for string in lls:
            new_string = string.replace("Â", "")
            if any(c in new_string for c in decomposition):
                continue
            range_matches = re.findall(
                r"([+-]?\d+[\.\d+]*)\s*-\s*([+-]?\d+[\.\d+]*)\s*°([CF])", new_string
            )
            if range_matches:
                range_matches = re.findall(
                    r"([+-]?\d+[\.\d+]*)\s*-\s*([+-]?\d+[\.\d+]*)\s*°([CF])", new_string
                )  # separates temps into ['100', 'C']
                if "C" in range_matches[0][-1]:
                    c_temp.append(range_matches[0][0])
                    c_temp.append(range_matches[0][1])

                elif "F" in range_matches[0][-1]:
                    f_temp.append(range_matches[0][0])
                    f_temp.append(range_matches[0][1])
                continue

            matches = re.findall(
                r"([+-]?\d*(\.\d+)*)\s*°([CF])", new_string
            )  # separates temps into ['100', 'C']
            if not matches:
                continue
            elif "C" in matches[0][-1]:
                c_temp.append(matches[0][0])
            elif "F" in matches[-1]:
                f_temp.append(matches[0][0])
        return c_temp, f_temp

    def fahrenheit_to_celsius(lls):
        """

        :param lls:
        :return temp_ls_in_c:
        """
        temp_ls_in_c = []
        for temp in lls:
            temp_in_c = (temp - 32) * 5 / 9
            temp_ls_in_c.append(temp_in_c)
        return temp_ls_in_c

    def minmax(lls):
        """This function is for testing purposes to check the range of temperature values supplied are within a certain
        limit of each other"""
        min_val = min(lls)
        max_val = max(lls)
        val_range = max_val - min_val
        if val_range > 10:
            print("range greater than 10")

        else:
            print("range lower than 10")

    celsius_temps_str, fahrenheit_temps_str = extract_temp(lls)

    celsius_temps = []
    for i in celsius_temps_str:
        try:
            celsius_temps.append(float(i))
        except ValueError:
            continue

    fahrenheit_temps = []
    for i in fahrenheit_temps_str:
        try:
            fahrenheit_temps.append(float(i))
        except ValueError:
            continue

    converted_temps = fahrenheit_to_celsius(fahrenheit_temps)

    all_temps = celsius_temps + converted_temps

    if all_temps:
        median_temp = median(all_temps)
        rounded_median_temp = round(median_temp, 1)
        return rounded_median_temp
    else:
        pass


def density_eval(lls):
    """
    :param lls:
    :return density:"""

    density_ls = []
    wrong_props = [
        "decomposes",
        "Decomposes",
        "decomposition",
        "Decomposition",
        "vapour",
        "Vapour",
    ]

    for string in lls:
        if any(c in string for c in wrong_props):
            continue
        new_string = string.replace("Â", "")
        temp_matches = re.findall(r"([-+]?\d*\.\d+|\d+)\s?Â*°([CF])", new_string)
        if not temp_matches:
            density_matches = re.findall(r"[+]?\d+\.\d+", string)
            if density_matches:
                density_float = float(density_matches[0])
                density_ls.append(density_float)
        elif "C" in temp_matches[0][-1]:
            temp = float(temp_matches[0][0])
            if 20 <= temp <= 25:
                density_matches = re.findall(r"[+]?\d+\.\d+", string)
                if density_matches:
                    density_float = float(density_matches[0])
                    density_ls.append(density_float)

        elif "F" in temp_matches[0][-1]:
            temp = float(temp_matches[0][0])
            if 68 <= temp <= 77:
                density_matches = re.findall(r"[+]?\d+\.\d+", string)
                if density_matches:
                    density_float = float(density_matches[0])
                    density_ls.append(density_float)

        density_ls = [c for c in density_ls if c < 15]
        if density_ls:
            density = median(density_ls)
            rounded_density = round(density, 3)
            return rounded_density


def compound_assertions(
    cas,
    cid,
    common_name,
    inchi,
    inchi_key,
    smiles_str,
    molec_formula,
    mw,
    hazard_codes,
    conc,
    density,
    bp,
    mp,
    flash_point,
    autoignition_temp,
):
    """Asserts the properties are as expected
    :param cas:
    :param cid:
    :param name:
    :param inchi:
    :param inchi_key:
    :param smiles_str:
    :param molec_formula:
    :param mw:
    :param hazard_codes:
    :param density:
    :param bp:
    :param mp:
    :param flash_point:
    :param autoignition_temp::
    :return:
    """

    assert 7 <= len(cas) <= 13, "invalid CAS length"
    required_entries = (cas, cid, common_name, inchi, molec_formula, mw)
    for entry in required_entries:
        if entry is None:
            print(f"{entry=}")
        # assert entry is not None, f"The required entry {entry} is set to None"

    assert mw >= 0, "Molecular weight is negative "

    temperatures = (bp, mp, flash_point, autoignition_temp)
    for temperature in temperatures:
        if temperature is not None:
            assert temperature > -273.15, "Temperature below absolute zero"

    if density is not None:
        assert (
            15 > density
        ), "Density higher than 15, expected highest is ~14 for Mercury"

    if hazard_codes != "No hazard codes found":
        split_codes = hazard_codes.split("-")
        for code in split_codes:
            if code:
                assert code[0] == "H"


def add_compound(
    cas,
    cid,
    common_name,
    inchi,
    inchi_key,
    smiles_str,
    molec_formula,
    mw,
    hazard_codes,
    conc,
    density,
    bp,
    mp,
    flash_point,
    autoignition_temp,
    model: Union[models.Compound, models.UpdateCompound] = models.Compound,
) -> models.Compound:
    """
    Creates a compound in the db session.
    :param cas:
    :param cid:
    :param name:
    :param inchi:
    :param inchi_key:
    :param smiles_str:
    :param molec_formula:
    :param mw:
    :param hazard_codes:
    :param density:
    :param bp:
    :param mp:
    :param flash_point:
    :param autoignition_temp::
    :param bind: the database to bind to.
    :return:
    """
    return model.create(
        commit=False,
        cas=cas,
        cid=cid,
        name=common_name,
        inchi=inchi,
        inchikey=inchi_key,
        smiles=smiles_str,
        molec_formula=molec_formula,
        molec_weight=mw,
        hphrase=hazard_codes,
        concentration=conc,
        density=density,
        boiling_point=bp,
        melting_point=mp,
        flash_point=flash_point,
        autoignition_temp=autoignition_temp,
    )


def extract_from_pubchem(
    compounds_file: str,
    limit: int,
    model: Union[models.Compound, models.UpdateCompound] = models.Compound,
):
    """
    The main function: This goes through the pubchem database 1 record at a time and if a CAS number is
    present then the compound will be added to the database

    Args:
        compounds_file: Path to the compounds file
        limit: Limit of compounds to process.
        bind: Database bind to input data to.

    :return:
    """
    print(
        f"Extracting compound data from pubchem download with a compound limit of: {limit}"
    )
    context = et.iterparse(
        compounds_file,
        events=("end",),
        tag="{http://pubchem.ncbi.nlm.nih.gov/pug_view}Record",
        encoding="ISO-8859-1",
    )
    for idx, (event, records) in enumerate(context):  # event is end of record
        cas = get_cas_num(records)
        inchi_key = get_identifier(records, "InChI Key")
        cas_eval = cas_verification(cas, model)
        if cas_eval is False:
            continue
        cid = get_cid(records)  # CID is the PubChem ID and hence the Unique identifier
        common_name = get_compound_name(records)
        inchi = get_identifier(records, "InChI")
        smiles_str = get_smiles(inchi)
        molec_formula = get_mol_formula(records)
        mw = get_mol_weight(records)
        density = get_phys_prop(records, "Density")
        bp = get_phys_prop(records, "Boiling Point")
        mp = get_phys_prop(records, "Melting Point")
        flash_point = get_phys_prop(records, "Flash Point")
        autoignition_temp = get_phys_prop(records, "Autoignition Temperature")
        echa_ref = find_echa(records)
        hazard_codes = get_hazards(records, echa_ref)
        conc = None

        if hazard_codes is None:
            hazard_codes = "No hazard codes found"

        # None values are suitable for float entries. '' is required instead of None for str entries for missing data.
        compound_assertions(
            cas,
            cid,
            common_name,
            inchi,
            inchi_key,
            smiles_str,
            molec_formula,
            mw,
            hazard_codes,
            conc,
            density,
            bp,
            mp,
            flash_point,
            autoignition_temp,
        )
        add_compound(
            cas,
            cid,
            common_name,
            inchi,
            inchi_key,
            smiles_str,
            molec_formula,
            mw,
            hazard_codes,
            conc,
            density,
            bp,
            mp,
            flash_point,
            autoignition_temp,
            model,
        )
        if idx > limit:
            break


def seed_solvent_data():
    """
    Seed solvent data to the database.
    """
    # add xylene to compound db
    models.Compound.create(commit=False, cid=-1, cas="1330-20-7", name="Xylenes")
    # get each row of solvent data from solvent csv
    with open(compound_database_dir / "solvents.csv") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")
        for idx, row in enumerate(csv_reader):
            if idx == 0:
                continue
            name = row[0]
            flag = row[1]
            hazard = row[2]
            cas = row[3]

            compound_db_match = (
                db.session.query(models.Compound)
                .filter(models.Compound.cas == cas)
                .first()
            )
            if compound_db_match is None:
                print(f"no match for: {name}")
            else:
                print(f"Match found between:{name} and {compound_db_match.name}")

            if not (
                name.isupper() or name in ["t-Butanol", "n-Butanol", "n-Butylacetate"]
            ):
                name = name.capitalize()
            if compound_db_match:
                models.Solvent.create(
                    commit=False,
                    name=name,
                    flag=flag,
                    hazard=hazard,
                    compound=[compound_db_match],
                    time_of_creation=current_time,
                )
            else:
                models.Solvent.create(
                    commit=False,
                    name=name.capitalize(),
                    flag=flag,
                    hazard=hazard,
                    time_of_creation=current_time,
                )


def seed_hazard_code_data():
    """
    Seeds hazard code data to the database.
    """
    haz_con = sql.connect(compound_database_dir / "hazard_codes.db")
    haz_con.row_factory = sql.Row
    haz_cur = haz_con.cursor()
    haz_cur.execute("SELECT * from hazard_codes_table")
    haz_rows = haz_cur.fetchall()
    haz_con.close()
    for hazard in haz_rows:
        code = hazard[0]
        phrase = hazard[1]
        category = hazard[2]
        models.HazardCode.create(
            commit=False, code=code, phrase=phrase, category=category
        )


def seed_element_sustainability_data():
    """
    Seeds element sustainability data to the database.
    """
    with open(compound_database_dir / "elements_sustainability.csv") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")
        for idx, row in enumerate(csv_reader):
            if idx == 0:
                continue
            element_name = row[0]
            chemical_symbol = row[1]
            remaining_supply = row[2]
            colour = row[3]
            models.Element.create(
                commit=False,
                name=element_name,
                symbol=chemical_symbol,
                remaining_supply=remaining_supply,
                colour=colour,
            )


def seed_role_types():
    """
    Seeds role types to the database.
    """
    models.Role.create(
        commit=False,
        name="Admin",
        role_description="The site administrator with access to the administration tabs",
    )
    models.Role.create(
        commit=False,
        name="Standard",
        role_description="The role assigned to new users signing up (PIs, SRs, SMs)",
    )


def seed_predefined_data():
    """
    Seeds predefined data to the database.
    """
    # Adds an institution then reads the yaml and adds the contents from there
    institution = models.Institution.create(
        commit=False, name="Test User Institution", time_of_creation=current_time
    )

    # Adding predefined data from the yaml file.
    seed_data_yaml_filepath = os.path.join(os.path.dirname(basedir), 'Webapp', 'seed_data.yaml')
    users_to_add = read_yaml(["predefined_users"], seed_data_yaml_filepath)
    try:
        workgroups_to_add = read_yaml(["predefined_workgroups"], seed_data_yaml_filepath)
    except Exception:
        workgroups_to_add = None
    for user in users_to_add.values():
        user_attributes = ["username", "email", "password", "fullname"]
        for attribute in user_attributes:
            if user[attribute] == "":
                print(
                    "\nWarning\nDatabase creation aborted\nYou must define admin account credentials in "
                    "Webapp/configs.yaml and rerun 'pubchem_download.py to setup the application correctly."
                )
                exit()
        p = models.Person.create(commit=False)

        role = "Standard"
        if user["admin"] is True:
            role = "Admin"
        role = db.session.query(models.Role).filter(models.Role.name == role).first()
        models.User.create(
            commit=False,
            username=user["username"],
            email=user["email"],
            fullname=user["fullname"],
            person=p.id,
            password_hash=models.User.set_password(user["password"]),
            role=role.id,
        )
    if workgroups_to_add:
        for workgroup in workgroups_to_add.values():
            workgroup_users = {}
            for user_type in [
                "principal_investigators",
                "senior_researchers",
                "standard_members",
            ]:
                workgroup_user_emails = workgroup[user_type].split(";")
                persons = (
                    db.session.query(models.Person)
                    .join(models.User)
                    .filter(models.User.email.in_(workgroup_user_emails))
                    .all()
                )
                workgroup_users[user_type] = persons
            workgroup_obj = models.WorkGroup.create(
                commit=False,
                name=workgroup["name"],
                institution=institution.id,
                principal_investigator=workgroup_users["principal_investigators"],
                senior_researcher=workgroup_users["senior_researchers"],
                standard_member=workgroup_users["standard_members"],
            )

            # add workbooks which belong to the workgroup
            workbooks_to_add = workgroup["workbooks"]
            for workbook in workbooks_to_add.values():
                workbook_members_emails = workbook["members"].split(";")
                persons = (
                    db.session.query(models.Person)
                    .join(models.User)
                    .filter(models.User.email.in_(workbook_members_emails))
                    .all()
                )
                models.WorkBook.create(
                    commit=False,
                    name=workbook["name"],
                    abbreviation=workbook["abbreviation"],
                    users=persons,
                    group=workgroup_obj.id,
                )
    db.session.commit()
