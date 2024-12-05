import os

from pony.orm import db_session, select
from pony.orm.core import TransactionIntegrityError
from sources import app, db
from datetime import datetime
import time
import pytz
import json
from utilities import read_yaml

DB_PROVIDER = read_yaml(['database_configurations', 'active_db'])

def test_database_create(db):
    if 'UNIT_TEST' not in os.environ:
        print("Add Environmental variable: 'UNIT_TEST=1' when running a unit test")
        exit()

    dropped_tables_ls = ['Reaction', 'User', 'Person', 'WorkBook', 'WorkGroup', 'Institution', 'NovelCompound',
                         'Reaction_Reaction', 'Person_WorkBook', 'Person_WorkGroup', 'Person_WorkGroup_2',
                         'Person_WorkGroup_3', 'WorkGroup_request', 'Notification', 'WBStatusRequest',
                         'WGStatusRequest', 'NewsItem', 'ReactionNote', 'ReactionDataFile']

    try:
        for table in dropped_tables_ls:
            db.drop_table(table, with_all_data=True)
    except:
        pass
    db.disconnect()
    db_config = app.config['PONY_DATABASE']
    db.provider = None
    db.schema = None
    db_config['provider'] = DB_PROVIDER
    # app.config.from_object('config.TestConfig')
    db.bind(db_config)
    #db.bind(provider=db_config['provider'], filename=db_config['filename'], create_db=db_config['create_db'])
    db.generate_mapping(create_tables=True)  # Makes tables based off the structure imported from models
    try:
        with db_session:
            p1 = db.Person()
            # Making a senior researcher
            p2 = db.Person()
            # Making a standard user
            p3 = db.Person()
            p4 = db.Person()
            p5 = db.Person()
            p6 = db.Person()

            u1 = db.User(username='PI_Test', email='PI@test.com', fullname='Pat Inglis',
                         person=p1, password_hash=db.User.set_password('PI_login'))
            u2 = db.User(username='SR_Test', email='SR@test.com', fullname='Sam Reed',
                         person=p2, password_hash=db.User.set_password('SR_login'))
            u3 = db.User(username='SM_Test', email='SM@test.com', fullname='Susan Matthews',
                         person=p3, password_hash=db.User.set_password('SM_login'))
            u4 = db.User(username='BB_Test', email='BB@test.com', fullname='Bob Brown',
                         person=p4, password_hash=db.User.set_password('BB_login'))
            u5 = db.User(username='SM1_Test', email='SM1@test.com', fullname='Stewart Mathers',
                         person=p5, password_hash=db.User.set_password('SM1_login'))
            u6 = (db.User(username='admin', email="admin@test.com", fullname="Admin Istrator",
                          person=p6, password_hash=db.User.set_password("admin_login"), role=db.Role[1]))
            institution1 = db.Institution(name='Test University')
            institution2 = db.Institution(name='University of Testing')
            institution3 = select(x for x in db.Institution if x.name == 'Test User Institution').first()
            if institution3 is None:
                institution3 = db.Institution(name="Test User Institution")
            # Making a PI

            # Making a workgroup owned by the SR and part of the institution
            wg = db.WorkGroup(name='Test-Workgroup', institution=institution1, principal_investigator=p1,
                              senior_researcher=p2, standard_member=[p3, p5])
            wg3 = db.WorkGroup(name='Test-Workgroup-2', institution=institution1, principal_investigator=p1,
                               senior_researcher=p2)
            wg4 = db.WorkGroup(name='Test-Workgroup-3', institution=institution2, principal_investigator=p4)

            # Making a workbook as part of the above workgroup
            wb = db.WorkBook(name='Test-Workbook', group=wg, users=[p1, p2, p3], abbreviation='TW1')
            wb3 = db.WorkBook(name='Test-Workbook2', group=wg, users=[p1, p2], abbreviation='TW2')
            wb4 = db.WorkBook(name='Test-Workbook3', group=wg3, users=[p1, p3], abbreviation='TW3')
            wb5 = db.WorkBook(name='Test-Workbook4', group=wg3, users=[p1, p2], abbreviation='TW4')
            wb6 = db.WorkBook(name='Test-Workbook5', group=wg4, users=p4, abbreviation='TW5')



    except TransactionIntegrityError:
        pass
    try:
        with db_session:
            wb3 = select(x for x in db.WorkBook if x.name == 'Test-Workbook2').first()
            wb6 = select(x for x in db.WorkBook if x.name == 'Test-Workbook5').first()
            p1 = select(x for x in db.Person if x.user.email == 'BB@test.com').first()
            # reaction1 = db.Reaction(creator=p1, time_of_creation='2021-10-06 14:51:43.311066',
            #                         name='nitro reduction', description='iron catalysed nitro reduction',
            #                         reactants='2778', products='1818', reagents=['9165', '9666'], reaction_smiles='[O-][N+](=O)C1=CC=CC=C1>>NC1=CC=CC=C1',
            #                         solvent='206', workbooks=wb, reaction_table_data=reaction_table, summary_table_data=summary_table)
            time_of_creation = datetime.now(pytz.timezone('Europe/London')).replace(tzinfo=None)
            # get compound ids
            nitrobenzene_id = select(x.id for x in db.Compound if x.name == 'Nitrobenzene').first()
            water_id = select(x.id for x in db.Compound if x.name == 'Water').first()
            methane_id = select(x.id for x in db.Compound if x.name == 'Methane').first()
            octane_id = select(x.id for x in db.Compound if x.name.lower() == 'octane').first()
            ethanol_id = select(x.id for x in db.Compound if x.name == 'Ethanol').first()
            thf_id = select(x.id for x in db.Compound if x.name == 'Tetrahydrofuran').first()
            aniline_id = select(x.id for x in db.Compound if x.name == 'Aniline').first()

            nitrobenzene_smiles = select(x.smiles for x in db.Compound if x.name == 'Nitrobenzene').first()
            water_smiles = select(x.smiles for x in db.Compound if x.name == 'Water').first()
            methane_smiles = select(x.smiles for x in db.Compound if x.name == 'Methane').first()
            octane_smiles = select(x.smiles for x in db.Compound if x.name.lower() == 'octane').first()
            ethanol_smiles = select(x.smiles for x in db.Compound if x.name == 'Ethanol').first()
            thf_smiles = select(x.smiles for x in db.Compound if x.name == 'Tetrahydrofuran').first()
            aniline_smiles = select(x.smiles for x in db.Compound if x.name == 'Aniline').first()

            # data for saving a reaction
            reaction_table = json.dumps({'reaction_smiles': '[O-][N+](=O)C1=CC=CC=C1>>NC1=CC=CC=C1',
                              'reaction_name': 'nitro reduction2',
                              'reaction_description': 'iron catalysed nitro reduction',
                              # reactant data
                              'reactant_ids': [nitrobenzene_id, water_id],
                              'reactant_masses': ['123', '180'],
                              'reactant_masses_raw': ['123', '180'],
                              'reactant_amounts': ['1.00', '10.0'],
                             'reactant_amounts_raw': ['1.00', '10.0'],
                              'reactant_volumes': ['0.10', '0.02'],
                              'reactant_densities': ['1.204', '0.995'],
                              'reactant_volumes_raw': ['0.10', '0.02'],
                              'reactant_equivalents': ['1', '10'],
                              'reactant_concentrations': ['-', '-'],
                              'reactant_physical_forms': ['2', '6'],
                              'limiting_reactant_table_number': 1,

                              'reagent_ids': [methane_id, octane_id],
                              'reagent_names': ['Methane', 'Octane'],
                              'reagent_amounts': ['0.20', '2.60'],
                              'reagent_amounts_raw': ['0.20', '2.60'],
                              'reagent_equivalents': ['0.2', '2.6'],
                              'reagent_molecular_weights': ['16.043', '114.23'],
                              'reagent_densities': ['', '0.703'],
                              'reagent_concentrations': ['', ''],
                              'reagent_masses': ['3.21', '297'],
                              'reagent_masses_raw': ['3.2054', '296.706'],
                              'reagent_volumes': ['', '0.42'],
                              'reagent_volumes_raw': ['', '0.422475'],
                              'reagent_physical_forms': ['1', '1'],
                              'reagent_hazards': ['H220', 'H225-H302-H304-H312-H315-H319-H336-H400-H410'],

                              'solvent_names': ['Ethanol', 'Tetrahydrofuran'],
                              'solvent_concentrations': ['0.50', '1.00'],
                              'solvent_ids': [ethanol_id, thf_id],
                              'solvent_volumes': ['2', '1'],
                              'solvent_physical_forms': ['2', '6'],
                              'solvent_hazards': ['H225-H302-H319-H371', 'H225-H302-H319-H335-H351'],

                              'product_ids': [aniline_id, water_id],
                              'product_masses': ['93.0', '18.0'],
                              'product_masses_raw': ['93.0', '18.0'],
                              'product_amounts': ['1.00', '1.00'],
                              'product_amounts_raw': ['1.00', '1.00'],
                              'product_mass_units': 'g',
                              'product_physical_forms': ['7', '6'],
                              'main_product_table_number': 7,
                              'main_product': 1,

                              'amount_units': 'mmol',
                              'mass_units': 'mg',
                              'volume_units': 'mL',
                              'solvent_volume_units': 'mL',
                              'product_amount_units': 'mmol'
                              })

            summary_table = json.dumps({'real_product_mass': '90',
                             'unreacted_reactant_mass': '2',
                             'reaction_temperature': '70',
                             'element_sustainability': '3',
                             'batch_flow': 'Batch',
                             'isolation_method': '2',
                             'catalyst_used': 'Catalyst or enzyme',
                             'catalyst_recovered': 'Not recovered catalyst',
                             'radio_buttons': ['hplc', 'massSpec'],
                             'custom_protocol1': '',
                             'custom_protocol2': '',
                             'other_hazards_text': 'Weigh out reactants in fumehood',
                             'researcher': '',
                             'supervisor': ''})


            reaction2 = db.Reaction(creator=p1, time_of_creation=time_of_creation, reaction_id='TW6-001',
                                    name='a nitro reduction2', description='iron catalysed nitro reduction',
                                    reactants=[nitrobenzene_smiles, water_smiles], products=[aniline_smiles, water_smiles],
                                    reagents=[methane_smiles, octane_smiles], solvent=[ethanol_smiles, thf_smiles], reaction_smiles='[O-][N+](=O)C1=CC=CC=C1.O>>NC1=CC=CC=C1.OCC',
                                    workbooks=wb6, reaction_table_data=reaction_table, summary_table_data=summary_table, complete="not complete", status="active")

            time.sleep(1)
            time_of_creation2 = datetime.now(pytz.timezone('Europe/London')).replace(tzinfo=None)
            reaction3 = db.Reaction(creator=p1, time_of_creation=time_of_creation2, reaction_id='TW6-002',
                                    name='nitro reduction2', description='iron catalysed nitro reduction',
                                    reactants=['2778', '255'], products=['1818', '255'],
                                    reagents=['9140', '9666'], solvent=['167', '206'],
                                    reaction_smiles='[O-][N+](=O)C1=CC=CC=C1.O>>NC1=CC=CC=C1.O',
                                    workbooks=wb6, reaction_table_data=reaction_table, summary_table_data=summary_table,
                                    complete="not complete", status="active")
            inchi_1 = 'InChI=1S/C21H15NO/c23-21(15-8-2-1-3-9-15)22-20-18-12-6-4-10-16(18)14-17-11-5-7-13-19(17)20/h1-14H,(H,22,23)'
            inchi_2 = 'InChI=1S/C10H16N4/c1-2-9(1)10-7-12-14(8-10)13-5-3-11-4-6-13/h7-9,11H,1-6H2'
            novel_reactant1 = db.NovelCompound(name='Starting material', cas='', workbook=wb6, molec_formula='C21H15NO',
                                               hphrase='Unknown', density=None, concentration=None, molec_weight=297.12,
                                               smiles='O=C(NC1=C2C=CC=CC2=CC2=C1C=CC=C2)C1=CC=CC=C1', InChI=inchi_1, InChIKey='MRGYEWYAVKQUEA-UHFFFAOYSA-N')
            novel_reactant2 = db.NovelCompound(name='Cyclopropane-series1', workbook=wb3, hphrase='Unknown', density=None,
                                               molec_formula='C10H16N4', concentration=None, molec_weight=192.14,
                                               smiles='C1CC1C1=CN(N=C1)N1CCNCC1', InChI=inchi_2, InChIKey='WBTXPQTZXVITMN-UHFFFAOYSA-N')
    except TransactionIntegrityError:
            pass

    return db


db = test_database_create(db)
