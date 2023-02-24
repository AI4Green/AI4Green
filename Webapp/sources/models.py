"""
This module contains database models that represent the data
that will be stored in the reaction database by a collection of classes
"""
from pony.orm import Optional, Set, StrArray, composite_key, Json, Required, PrimaryKey
from datetime import datetime
from flask_login import UserMixin
from werkzeug.security import generate_password_hash, check_password_hash
from sources import login


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
        role = Required('Role', default=2)  # 2 is a normal user, 1 is admin
        person = Required('Person')
        hazard_colors = Required(Json, default={"Recommended": "#00ff00", "Problematic": "#ffff00",
                                                "Hazardous": "#ff0000", "HighlyHazardous": "#8B0000",
                                                "Recommended_text": "#000000", "Problematic_text": "#000000",
                                                "Hazardous_text": "#000000", "HighlyHazardous_text": "#ffffff"
                                                })

        """Password hashing is implemented 
		by the two following methods"""

        def check_password(self, password):
            return check_password_hash(self.password_hash, password)

        @classmethod
        def set_password(cls, password):
            return generate_password_hash(password)

        """"""

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

    @login.user_loader
    def load_user(user_id):
        """Flask-Login keeps track of the logged in user
		by storing its unique identifier in Flask's user
		session. Each time the logged-in user navigates
		to a new page, Flask-Login retrieves the ID of
		the user from the session, and then loads that
		user into memory by the user loader function"""
        return db.User.get(id=user_id)
