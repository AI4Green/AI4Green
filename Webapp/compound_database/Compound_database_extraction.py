"""
The code here is for constructing the chemical database using the XML from PubChem
"""
from datetime import datetime
from lxml import etree as et
from statistics import mode, median
from rdkit import Chem
import os
from flask_login import UserMixin
from pony.orm import Database, Required, Optional, StrArray, select, db_session, Set, composite_key, Json, PrimaryKey
import re
from chemspipy import ChemSpider
import sqlite3 as sql
import csv
# import cirpy
from utilities import read_yaml
from werkzeug.security import generate_password_hash, check_password_hash

ns = {'xmlns': 'http://pubchem.ncbi.nlm.nih.gov/pug_view',  # Defines the url used in the xml tags
      'xs': 'http://pubchem.ncbi.nlm.nih.gov/pug_view'}


def define_database(db):
    """This database defines tables for: Users, Compounds, Solvents, and linked tables that were previously in the reaction
    database: Person, Reaction, Novel Compound, Institution, Workbook, Workgroup"""

    class User(db.Entity, UserMixin):
        """The User class inherits from db.entity which is used as a base class
		 for all entities stored in the database and UserMixin which provides
		 default implementations for all of these properties and methods
		 This class defines attributes as class variables where the kind
		 and type of attribute and attribute options are defined."""
        _table_ = "User"
        username = Required(str, unique=True)
        email = Required(str, unique=True)
        fullname = Required(str)
        password_hash = Required(str)
        role = Required('Role', default=2)
        person = Required('Person')
        hazard_colors = Required(Json, default={"Recommended": "#00ff00", "Problematic": "#ffff00",
                                                "Hazardous": "#ff0000", "HighlyHazardous": "#8B0000",
                                                "Recommended_text": "#000000", "Problematic_text": "#000000",
                                                "Hazardous_text": "#000000", "HighlyHazardous_text": "#ffffff"
                                                })

        def check_password(self, password):
            return check_password_hash(self.password_hash, password)

        @classmethod
        def set_password(cls, password):
            return generate_password_hash(password)

    class Role(db.Entity):
        _table_ = "Role"
        name = Required(str)
        role_description = Required(str)
        users = Set(User)

    class Solvent(db.Entity):
        _table_ = "Solvent"
        name = PrimaryKey(str)
        flag = Required(int, nullable=True)
        hazard = Optional(str, nullable=True)
        compound = Optional('Compound', reverse='solvent')
        novel_compound = Optional('NovelCompound', reverse='solvent')


    class Compound(db.Entity):
        _table_ = "Compound"
        CID = Required(int, unique=True, nullable=True)
        cas = Required(str, unique=True)
        name = Required(str)
        smiles = Optional(str, nullable=True)
        InChI = Optional(str, nullable=True)
        InChIKey = Optional(str, nullable=True)

        molec_formula = Optional(str, nullable=True)
        density = Optional(float, nullable=True)
        concentration = Optional(float, nullable=True)

        boiling_point = Optional(float, nullable=True)  # bp = boiling point
        melting_point = Optional(float, nullable=True)  # mp = melting point
        flash_point = Optional(float, nullable=True)
        autoignition_temp = Optional(float, nullable=True)
        molec_weight = Optional(float, nullable=True)  # mw = molecular weight
        state = Optional(str, nullable=True)
        form = Optional(str, nullable=True)
        hphrase = Optional(str, nullable=True)  # in future this should be an optional StrArray

        safety_score = Optional(float, nullable=True)
        health_score = Optional(float, nullable=True)
        enviro_score = Optional(float, nullable=True)
        econom_score = Optional(float, nullable=True)
        error_report = Set("CompoundDataErrorReport", nullable=True)
        solvent = Optional('Solvent')



    class HazardCode(db.Entity):
        _table_ = "HazardCode"
        code = Required(str)
        phrase = Required(str)
        category = Optional(str, nullable=True)

    class Reaction(db.Entity):
        """
        A Reaction is defined by;
        name: identifier
        date: date of creation of the reaction
        date_reaction: date when the reaction was performed
        precursor_reaction: links to the reaction in the database that preceded it
        successor_reaction: links to the reactions that follow on from this one
        green_metric: the score given for this reaction
        reactants: links to a set of reactants for this reaction
        products: links to a set of products for this reaction
        reagents: links to a set of reagents for this reaction
        solvent: links to a set of solvents for this reaction
        name_book_comp: this is the Unique PrimaryKey id for the reaction based on the
            name of the reaction and the workbook it is performed in
        """
        _table_ = "Reaction"
        reaction_id = Required(str)
        name = Required(str)
        description = Optional(str)
        reaction_class = Optional(str)  # Not used yet
        creator = Required('Person')
        time_of_creation = Required(datetime, 6)
        time_of_update = Optional(datetime, 6)
        date_reaction = Optional(datetime, 6)
        precursor_reaction = Set('Reaction', table="Reaction_Reaction")
        successor_reaction = Set('Reaction', table="Reaction_Reaction")
        green_metric = Optional(float)
        workbooks = Required('WorkBook')
        reactants = Optional(StrArray)  # Required(StrArray)
        products = Optional(StrArray)  # Required(StrArray)
        reagents = Optional(StrArray)
        solvent = Optional(StrArray)
        reaction_table_data = Optional(Json)
        summary_table_data = Optional(Json)
        reaction_smiles = Optional(str)
        complete = Required(str)
        status = Required(str)
        name_book_comp = composite_key(workbooks, name)
        reaction_id_book_comp = composite_key(workbooks, reaction_id)

    class Person(db.Entity):
        """
		A Person is a user, and can have some attributes that give them further access
		rights within the Reaction Database

		email: Required unique identifier of the user
		fullname: Required
		workgroups: Optional Set of links to the workgroups (labs) that this user manages
		workgroup_principal_investigator: This is a set of links to the workgroups which this user is the PI/owner of
		workgroup_senior_researcher: This is a set of links to the workgroups which this user is a senior researcher in
		workgroup_standard_member: This is a set of links to the workgroups which this user is a standard member in
		workbook_user: This is a Set of links to the workbooks this user has access to
		"""
        _table_ = "Person"
        user = Optional(User)
        workgroup_principal_investigator = Set('WorkGroup', reverse='principal_investigator',
                                               table="Person_WorkGroup")  # PI is also the owner
        workgroup_principal_investigator_request = Set("WorkGroup_request")  # PI is also the owner
        person_notification = Set("Notification")
        person_status_change = Set("WGStatusRequest", reverse='person')
        person_status_change_pi = Set("WGStatusRequest", reverse='principal_investigator')
        person_status_change_wb = Set("WBStatusRequest", reverse='person')
        person_status_change_wb_sr_pi = Set("WBStatusRequest", reverse='pi_sr')
        workgroup_senior_researcher = Set('WorkGroup', reverse='senior_researcher', table="Person_WorkGroup_2")
        workgroup_standard_member = Set('WorkGroup', reverse='standard_member', table="Person_WorkGroup_3")
        workbook_user = Set('WorkBook', table="Person_WorkBook")
        reaction = Set('Reaction')

    class WorkGroup(db.Entity):
        """
        A workgroup is a Person's lab - if they run one.
        A workgroup thus has a name and person associated with it, that gives the Primary key
        A workgroup can have other owners - this makes sense for larger organisations
        workgroup_members is the list of Person's who are part of this workgroup
        workbooks are a set of projects (lab books) that the lab maintains
        """
        _table_ = "WorkGroup"
        name = Required(str)
        principal_investigator = Set(Person, table="Person_WorkGroup")  # The PI is the workgroup owner.
        senior_researcher = Set(Person, table="Person_WorkGroup_2")
        standard_member = Set(Person, table="Person_WorkGroup_3")
        WorkGroup_Status_Request = Set("WGStatusRequest")
        workbooks = Set('WorkBook')
        institution = Required('Institution')
        institution_group_name_comp = composite_key(name, institution)  # Must be unique
        approved = Optional(bool)
        request = Optional('WorkGroup_request', cascade_delete=True)

    class WorkGroup_request(db.Entity):
        """
        A workgroup is a Person's lab - if they run one.
        A workgroup thus has a name and person associated with it, that gives the Primary key
        A workgroup can have other owners - this makes sense for larger organisations
        workgroup_members is the list of Person's who are part of this workgroup
        workbooks are a set of projects (lab books) that the lab maintains
        """
        _table_ = "WorkGroup_request"
        name = Required(str)
        principal_investigator = Required(Person)  # The PI is the workgroup owner.
        info = Required(str)
        time = Required(datetime)
        status = Required(str)
        institution = Required('Institution')
        institution_group_name_comp = composite_key(name, institution)  # Must be unique
        workgroup = Required(WorkGroup)

    class WorkBook(db.Entity):
        """
        name is the name of the workbook
        group is the workgroup the workbook is associated
        users is a list of Persons that are allowed to access a workbook
        reactions is a list of reactions within the workbook
        the name and group is used to create the PrimaryKey
        """
        _table_ = "WorkBook"
        name = Required(str)
        abbreviation = Required(str)
        group = Required(WorkGroup)
        users = Set(Person, table="Person_WorkBook")
        reactions = Set(Reaction)
        novel_compound = Set('NovelCompound')
        WorkBook_Status_Request = Set("WBStatusRequest")
        name_group_comp = composite_key(name, group)  # Must be unique

    class NovelCompound(db.Entity):
        _table_ = "NovelCompound"
        name = Required(str)
        cas = Optional(str, nullable=True)
        smiles = Optional(str)
        InChI = Optional(str)
        InChIKey = Optional(str, nullable=True)

        molec_weight = Optional(float, nullable=True)
        molec_formula = Optional(str)
        density = Optional(float, nullable=True)
        concentration = Optional(float, nullable=True)

        boiling_point = Optional(float, nullable=True)
        melting_point = Optional(float, nullable=True)
        flash_point = Optional(float, nullable=True)
        autoignition_temp = Optional(float, nullable=True)
        state = Optional(str, nullable=True)
        form = Optional(str, nullable=True)
        hphrase = Optional(str, nullable=True)  # in future this should be a Optional(StrArray)

        safety_score = Optional(float)
        health_score = Optional(float)
        enviro_score = Optional(float)
        econom_score = Optional(float)

        workbook = Required(WorkBook)
        name_book_comp = composite_key(name, workbook)  # Compound name must be unique within workbook
        structure_book_comp = composite_key(InChI, workbook)  # Compound InChI must be unique within workbook
        solvent = Optional(Solvent, cascade_delete=True)

    class Institution(db.Entity):
        _table_ = "Institution"
        name = Required(str, unique=True)
        workgroups = Set(WorkGroup)
        workgroups_request = Set(WorkGroup_request)

    class CompoundDataErrorReport(db.Entity):
        _table_ = "CompoundDataErrorReport"
        compound_name = Required(str)
        compound = Required(Compound)
        error_type = Required(str)
        additional_info = Optional(str)
        # user = Required(User)
        time = Required(datetime)

    class Notification(db.Entity):
        _table_ = "Notification"
        person = Required(Person)
        type = Required(str)
        info = Required(str)
        time = Required(datetime)
        status = Required(str)
        WG = Optional(str)
        WB = Optional(str)
        wg_status_notification = Set("WGStatusRequest")
        wb_status_notification = Set("WBStatusRequest")
        wg_request = Optional(str)

    class WGStatusRequest(db.Entity):
        _table_ = "WGStatusRequest"
        principal_investigator = Required(Person)
        person = Required(Person)
        WG = Required(WorkGroup)
        current_role = Required(str)
        new_role = Required(str)
        time = Required(datetime)
        notification = Required(Notification)
        status = Required(str)

    class WBStatusRequest(db.Entity):
        _table_ = "WBStatusRequest"
        pi_sr = Required(Person)
        person = Required(Person)
        WB = Required(WorkBook)
        current_role = Required(str)
        new_role = Required(str)
        time = Required(datetime)
        notification = Required(Notification)
        status = Required(str)

    class Element(db.Entity):
        _table_ = "Element"
        name = Required(str)
        symbol = Required(str)
        remaining_supply = Required(str)
        colour = Required(str)

    class NewsItem(db.Entity):
        _table_ = "NewsItem"
        title = Required(str)
        message = Required(str)
        time = Required(datetime)


def open_database_compound(db_config):
    print(f"Making database with config:{db_config}")
    db = Database()
    define_database(db)
    db.bind(db_config)
    db.generate_mapping(create_tables=True)
    return db


def get_cas_num(root):
    """Returns the CAS number when given the root for a single record.
    :param root:
    :return cas_number, cas_number_ls:
    """

    for sections in root.findall('xmlns:Section/xmlns:Section', ns):
        # Finds where the CAS number is stored in the XML
        for headings in sections.findall('xmlns:TOCHeading', ns):
            if headings.text == 'CAS':
                # Saves the CAS numbers as a list
                cas_number_ls = \
                    [casNum.text for casNum in sections.findall('xmlns:Information/xmlns'
                                                                ':Value/xmlns:StringWithMarkup/xmlns:String',
                                                                ns)]
                # Finds the mode of the cas numbers.
                cas_number = mode(cas_number_ls)
                return cas_number
            else:
                return None  # returning both as none means the record will be skipped if no cas is found


def cas_verification(cas_number, db):
    """The purpose of this code is to determine if a CAS number is valid or not.
    If the CAS number is invalid or already present in the database the current compound record will not be added to
    the database and will skip to the next compound record by returning
    False.
    :param db:
    :param cas_number:
    :return:

    """
    with db_session:
        duplicate_cas_check = select(c for c in db.Compound if c.cas == cas_number).first()
        if duplicate_cas_check is not None:
            # print('duplicate cas - compound omitted from database')
            return False

    if cas_number is None:
        return False
    elif not 7 <= len(cas_number) <= 13:
        print('wrong length cas - compound omitted from database')
        return False
    else:
        return True


def get_compound_name(root):
    """This function retrieves the 'main' name for a compound that is in the record title of the xml file
    :param root:
    :return compound_name:

    """

    for tags in root.findall('xmlns:RecordTitle', ns):
        compound_name = tags.text
        return compound_name


def get_cid(root):
    """This function retrieves the pub chem compound identifier (CID) which is unique to a compound.
    Type should be either RecordNumber or RecordTitle
    :param root:
    :return record_info:
    """
    for tags in root.findall('xmlns:RecordNumber', ns):
        record_info = tags.text
        return record_info


def get_identifier(root, id_type):
    """This function is used to retrieve the InChI and InChI key from the xml file
    :param root:
    :param id_type:
    :return id:
    """
    for sections in root.findall('xmlns:Section/xmlns:Section', ns):

        for headings in sections.findall('xmlns:TOCHeading', ns):
            if headings.text == id_type:

                for id_values in sections.findall('xmlns:Information/xmlns:Value/xmlns:StringWithMarkup/xmlns:String',
                                                  ns):
                    ident = id_values.text
                    return ident


def get_smiles(inchi_str):
    """
    Function for converting the InChI string into a Smiles string
    :param inchi_str:
    :return:
    """

    try:
        mol = Chem.inchi.MolFromInchi(inchi_str, sanitize=True, removeHs=True, logLevel=None,
                                      treatWarningAsError=False)
        smiles_str = Chem.MolToSmiles(mol)
        return smiles_str
    except:
        pass

    try:
        mol = Chem.inchi.MolFromInchi(inchi_str, sanitize=False, logLevel=None,
                                      treatWarningAsError=False)
        smiles_str = Chem.MolToSmiles(mol)
        return smiles_str
    except:
        pass

    try:
        cs = ChemSpider('584QozVSGkkJn8YulyGGSthFGIgEB6PJ')
        smiles_str = (cs.convert(input=inchi_str, input_format='InChI', output_format='SMILES'))
        return smiles_str
    except:
        return 'smiles error'


def get_mol_weight(root):
    """This function returns the molecular weight when given the root for a single record
    :param root:
    :return molecular_weight:
    """

    for sections in root.findall('xmlns:Section', ns):
        # Finds where the molecular weight is stored in the xml
        for headings in sections.findall('xmlns:TOCHeading', ns):
            if headings.text == 'Molecular Weight':
                # Saves the molecular weight and converts it to a string.
                # molecular_weight_ls = [mol_weight.text for mol_weight in
                #                        sections.findall('xmlns:Information/xmlns:Value/xmlns:Number',
                #                                         ns)]
                molecular_weight_ls = [mol_weight.text for mol_weight in
                                       sections.findall('xmlns:Information/xmlns:Value/xmlns:StringWithMarkup/'
                                                        'xmlns:String',
                                                        ns)]
                molecular_weight = molecular_weight_ls[0]
                return float(molecular_weight)


def find_echa(root):
    """This function finds the ECHA reference number used for the hazard data.
    This means all hazard data extracted is the data supplied to PubChem by the ECHA
    :param root:
    :return ref_num:
    """

    for sections in root.findall('xmlns:Section/xmlns:Information', ns):

        for references in sections.findall('xmlns:Name', ns):
            if references.text == 'ECHA C&L Notifications Summary':
                ref_num_ls = [ref_nums.text for ref_nums in sections.findall('xmlns:ReferenceNumber', ns)]
                ref_num = ref_num_ls[0]
                return ref_num
        # alternative name
        for references in sections.findall('xmlns:Name', ns):
            if references.text == 'Signal':
                ref_num_ls = [ref_nums.text for ref_nums in sections.findall('xmlns:ReferenceNumber', ns)]
                ref_num = ref_num_ls[0]
                return ref_num


def get_hazards(root, echa_ref):
    """Returns hazard codes when given the root for a single record and the reference number for the echa data
    :param root:
    :param echa_ref:
    :return hazard_codes:
    """

    for sections in root.findall('xmlns:Section/xmlns:Information', ns):
        # Finds where the hazard codes are stored
        for headings in sections.findall('xmlns:Name', ns):
            if headings.text == 'GHS Hazard Statements':

                # Finds the ECHA supplied data using the reference number and saves them as a list
                for hazard_suppliers in sections.findall('xmlns:ReferenceNumber', ns):
                    if hazard_suppliers.text == echa_ref:

                        hazard_codes_ls = [haz_codes.text for haz_codes in
                                           sections.findall('xmlns:Value/xmlns:StringWithMarkup/xmlns:String',
                                                            ns)]

                        # Formats the hazard codes by separating the different codes then only saving the H codes.
                        hazard_codes = ""

                        for hazards in hazard_codes_ls:
                            hazards_text = hazards.split(" ")
                            ind_hazard_codes = hazards_text[0]
                            ind_hazard_codes = ind_hazard_codes.split(':')[0]
                            ind_hazard_codes += "-"
                            hazard_codes += ind_hazard_codes
                        hazard_codes = hazard_codes[:-1]
                        return hazard_codes


def get_mol_formula(root):
    """
    Extracts the molecular formula for a single record
    :param root:
    :return mol_formula:
    """

    for sections in root.findall('xmlns:Section', ns):

        for headings in sections.findall('xmlns:TOCHeading', ns):
            if headings.text == 'Molecular Formula':

                for values in sections.findall('xmlns:Information/xmlns:Value/xmlns:StringWithMarkup/xmlns:String',
                                               ns):
                    mol_formula = values.text
                    return mol_formula


def get_phys_prop(root, property_input):
    """This function is used to retrieve the physical properties from the xml file (temperatures, form and density)
    :param root:
    :param property_input:
    :return property:
    """

    for sections in root.findall('xmlns:Section/xmlns:Section', ns):

        for headings in sections.findall('xmlns:TOCHeading', ns):
            if headings.text == property_input:

                for entries in sections.findall('xmlns:Information/xmlns:Description',
                                                ns):
                    if entries.text == 'PEER REVIEWED':

                        phys_properties_ls = [values.text for values in sections.findall('xmlns:Information/xmlns'
                                                                                         ':Value/xmlns'
                                                                                         ':StringWithMarkup/xmlns'
                                                                                         ':String',
                                                                                         ns)]

                        temperatures = {"Boiling Point", "Melting Point", "Flash Point", "Autoignition Temperature"}

                        if phys_properties_ls:
                            if property_input in temperatures:
                                phys_property = temp_eval(phys_properties_ls)
                                return phys_property

                            elif property_input == 'Density':
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
        decomposition = ['decomposes', 'Decomposes', 'decomposition', 'Decomposition']
        for string in lls:
            new_string = string.replace("Â", "")
            if any(c in new_string for c in decomposition):
                continue
            range_matches = re.findall("([+-]?\d+[\.\d+]*)\s*-\s*([+-]?\d+[\.\d+]*)\s*°([CF])", new_string)
            if range_matches:
                range_matches = re.findall("([+-]?\d+[\.\d+]*)\s*-\s*([+-]?\d+[\.\d+]*)\s*°([CF])",
                                           new_string)  # separates temps into ['100', 'C']
                if 'C' in range_matches[0][-1]:
                    c_temp.append(range_matches[0][0])
                    c_temp.append(range_matches[0][1])

                elif 'F' in range_matches[0][-1]:
                    f_temp.append(range_matches[0][0])
                    f_temp.append(range_matches[0][1])
                continue

            matches = re.findall("([+-]?\d*(\.\d+)*)\s*°([CF])", new_string)  # separates temps into ['100', 'C']
            if not matches:
                continue
            elif 'C' in matches[0][-1]:
                c_temp.append(matches[0][0])
            elif 'F' in matches[-1]:
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

    conv_temps = fahrenheit_to_celsius(fahrenheit_temps)

    all_temps = celsius_temps + conv_temps

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
    wrong_props = ['decomposes', 'Decomposes', 'decomposition', 'Decomposition', 'vapour', 'Vapour']

    for string in lls:
        if any(c in string for c in wrong_props):
            continue
        new_string = string.replace("Â", "")
        temp_matches = re.findall("([-+]?\d*\.\d+|\d+)\s?Â*°([CF])", new_string)
        if not temp_matches:
            density_matches = re.findall("[+]?\d+\.\d+", string)
            if density_matches:
                density_float = float(density_matches[0])
                density_ls.append(density_float)
        elif 'C' in temp_matches[0][-1]:
            temp = float(temp_matches[0][0])
            if 20 <= temp <= 25:
                density_matches = re.findall("[+]?\d+\.\d+", string)
                if density_matches:
                    density_float = float(density_matches[0])
                    density_ls.append(density_float)

        elif 'F' in temp_matches[0][-1]:
            temp = float(temp_matches[0][0])
            if 68 <= temp <= 77:
                density_matches = re.findall("[+]?\d+\.\d+", string)
                if density_matches:
                    density_float = float(density_matches[0])
                    density_ls.append(density_float)

        density_ls = [c for c in density_ls if c < 15]
        if density_ls:
            density = median(density_ls)
            rounded_density = round(density, 3)
            return rounded_density


def compound_assertions(cas, cid, common_name, inchi, inchi_key, smiles_str, molec_formula, mw,
                        hazard_codes, conc,
                        density, bp, mp, flash_point, autoignition_temp):
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
        assert 15 > density, "Density higher than 15, expected highest is ~14 for Mercury"

    if hazard_codes != 'No hazard codes found':
        split_codes = hazard_codes.split('-')
        for code in split_codes:
            if code:
                assert code[0] == 'H'


@db_session
def add_compound(db, cas, cid, common_name, inchi, inchi_key, smiles_str, molec_formula, mw,
                 hazard_codes, conc,
                 density, bp, mp, flash_point, autoignition_temp):
    """Adds the compound to the db

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
    db.Compound(cas=cas, CID=cid, name=common_name, InChI=inchi, InChIKey=inchi_key,
                smiles=smiles_str, molec_formula=molec_formula, molec_weight=mw, hphrase=hazard_codes,
                concentration=conc,
                density=density,
                boiling_point=bp, melting_point=mp, flash_point=flash_point,
                autoignition_temp=autoignition_temp)


@db_session
def extract_from_pubchem(Compounds_file, db, limit):
    """
    The main function: This goes through the pubchem database 1 record at a time and if a CAS number is
    present then the compound will be added to the database

    :return:
    """
    print(f"Extracting compound data from pubchem download with a compound limit of: {limit}")
    context = et.iterparse(Compounds_file, events=("end",),
                           tag='{http://pubchem.ncbi.nlm.nih.gov/pug_view}Record', encoding='ISO-8859-1')
    for idx, (event, records) in enumerate(context):  # event is end of record

        cas = get_cas_num(records)
        inchi_key = get_identifier(records, 'InChI Key')
        cas_eval = cas_verification(cas, db)
        if cas_eval is False:
            continue
        cid = get_cid(records)  # CID is the PubChem ID and hence the Unique identifier
        common_name = get_compound_name(records)
        inchi = get_identifier(records, 'InChI')
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
            hazard_codes = 'No hazard codes found'

        # None values are suitable for float entries. '' is required instead of None for str entries for missing data.
        compound_assertions(cas, cid, common_name, inchi, inchi_key, smiles_str, molec_formula,
                            mw, hazard_codes, conc, density,
                            bp, mp, flash_point, autoignition_temp)
        add_compound(db, cas, cid, common_name, inchi, inchi_key, smiles_str, molec_formula,
                     mw, hazard_codes, conc, density,
                     bp, mp, flash_point, autoignition_temp)
        if idx > limit:
            break


@db_session
def migrate_solvent_data(db):
    print("Adding solvents")
    # add xylene to compound db
    db.Compound(CID=-1, cas='1330-20-7', name='Xylenes')
    # get each row of solvent data from solvent csv
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(BASE_DIR, 'solvents.csv')) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for idx, row in enumerate(csv_reader):
            if idx == 0:
                continue
            name = row[0]
            flag = row[1]
            hazard = row[2]
            cas = row[3]
            compound_db_match = select(x for x in db.Compound if x.cas == cas).first()
            if compound_db_match is None:
                print(f'no match for: {name}')
            else:
                print(f'Match found between:{name} and {compound_db_match.name}')
            if name.isupper() or name in ['t-Butanol', 'n-Butanol', 'n-Butylacetate']:
                db.Solvent(name=name, flag=flag, hazard=hazard, compound=compound_db_match)
            else:
                db.Solvent(name=name.capitalize(), flag=flag, hazard=hazard, compound=compound_db_match)



@db_session
def migrate_hazard_code_data(db):
    print("Adding hazard codes")
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    haz_con = sql.connect(os.path.join(BASE_DIR, "hazard_codes.db"))
    haz_con.row_factory = sql.Row
    haz_cur = haz_con.cursor()
    haz_cur.execute("SELECT * from hazard_codes_table")
    haz_rows = haz_cur.fetchall()
    haz_con.close()
    for hazard in haz_rows:
        code = hazard[0]
        phrase = hazard[1]
        category = hazard[2]
        # print(code, phrase, category)
        db.HazardCode(code=code, phrase=phrase, category=category)


@db_session
def migrate_element_sustainability_data(db):
    print("Adding element sustainability")
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(BASE_DIR, 'elements_sustainability.csv')) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for idx, row in enumerate(csv_reader):
            if idx == 0:
                continue
            element_name = row[0]
            chemical_symbol = row[1]
            remaining_supply = row[2]
            colour = row[3]
            db.Element(name=element_name, symbol=chemical_symbol, remaining_supply=remaining_supply, colour=colour)


@db_session
def add_role_types(db):
    print("Adding user roles")
    db.Role(name='Admin', role_description='The site administrator with access to the administration tabs')
    db.Role(name='Standard', role_description='The role assigned to new users signing up (PIs, SRs, SMs)')


@db_session
def add_predefined_data(db):
    # Adds an institution then reads the yaml and adds the contents from there
    institution1 = db.Institution(name='Test User Institution')
    # Adding predefined data from the yaml file.
    users_to_add = read_yaml(['predefined_users'])
    try:
        workgroups_to_add = read_yaml(['predefined_workgroups'])
    except:
        workgroups_to_add = None
    for user in users_to_add.values():
        user_attributes = ["username", "email", "password", "fullname"]
        for attribute in user_attributes:
            if user[attribute] == "":
                print("Please define admin account credentials in Webapp/configs.yaml and rerun 'pubchem_download.py "
                      "to setup the application correctly'")
                exit()
        p = db.Person()
        role = 'Standard'
        if user["admin"] is True:
            role = 'Admin'
        role = select(x for x in db.Role if x.name == role).first()
        db.User(username=user["username"], email=user["email"], fullname=user["fullname"],
                person=p, password_hash=db.User.set_password(user["password"]), role=role)
    if workgroups_to_add:
        for workgroup in workgroups_to_add.values():
            workgroup_users = {}
            for user_type in ["principal_investigators", "senior_researchers", "standard_members"]:
                workgroup_user_emails = workgroup[user_type].split(";")
                workgroup_user_ids = select(x.id for x in db.Person if x.user.email in workgroup_user_emails)[:]
                # for user_id in workgroup_user_ids:
                workgroup_users[user_type] = [db.Person[x] for x in workgroup_user_ids if db.Person[x]]
            workgroup_obj = db.WorkGroup(name=workgroup["name"], institution=institution1,
                                         principal_investigator=workgroup_users["principal_investigators"],
                                         senior_researcher=workgroup_users["senior_researchers"],
                                         standard_member=workgroup_users["standard_members"]
                                         )
            # add workbooks which belong to the workgroup
            workbooks_to_add = workgroup["workbooks"]
            for workbook in workbooks_to_add.values():
                workbook_members_emails = workbook["members"].split(";")
                workbook_members_ids = select(x.id for x in db.Person if x.user.email in workbook_members_emails)
                workbook_members = [db.Person[x] for x in workbook_members_ids if db.Person[x]]
                db.WorkBook(name=workbook['name'], abbreviation=workbook["abbreviation"], users=workbook_members,
                            group=workgroup_obj)
