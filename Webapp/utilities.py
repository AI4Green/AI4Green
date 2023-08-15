import os
import yaml
from shutil import copyfile
from typing import Optional

basedir = os.path.abspath(os.path.dirname(__file__))  # path to the user database
configs_yaml_filepath = os.path.join(os.path.dirname(basedir), 'Webapp', 'configs.yaml')


def read_yaml(target_keys: list, filepath: Optional[str] = configs_yaml_filepath):
    """
    target_keys: iterates through keys to find nested data in yaml config
    """
    with open(filepath, "r") as stream:
        try:
            target_data = yaml.safe_load(stream)
            for key in target_keys:
                target_data = target_data[key]
        except yaml.YAMLError as exc:
            print(exc)
        return target_data


def read_yaml_db_config(db_name):
    """Reads the yaml file and returns the db_config - applies changes if running a UNIT_TEST and sqlite filenames"""
    db_config = read_yaml(["database_configurations", db_name])
    if 'UNIT_TEST' in os.environ:
        change_to_unit_test_db(db_config)
    if db_config['provider'] == 'sqlite':
        # assign absolute path to filename
        db_config['filename'] = os.environ.get('DATABASE_URL') or os.path.join(basedir, 'test_database.sqlite')
    return db_config


def change_to_unit_test_db(db_config):
    """Selects appropriate function to edit db_config depending whether sqlite or postgres database is active_db"""
    if db_config['provider'] == 'sqlite':
        make_unit_test_sqlite_db(db_config)
    elif db_config['provider'] == 'postgres':
        use_unit_test_postgres_db(db_config)


def make_unit_test_sqlite_db(db_config):
    """Makes a copy of the sqlite database and uses this for unit tests"""
    copyfile(os.path.join(basedir, 'test_database.sqlite'), os.path.join(basedir, 'unit_test_database.sqlite'))
    db_config['filename'] = os.path.join(basedir, 'unit_test_database.sqlite')


def use_unit_test_postgres_db(db_config):
    """Remakes an instance of the database for unit tests if required."""
    from compound_database.pubchem_download import create_database, check_if_existing_postgres_db, \
        connect_to_postgres_db
    # change to connect to the unit test database
    db_config['database'] = 'unit_test'
    # check if already present
    cursor = connect_to_postgres_db(db_config)
    db_present = check_if_existing_postgres_db(cursor, db_config)
    # if not present then create
    if not db_present:
        create_database(db_config)


def db_config_for_print(db_config):
    """Prevents password printing to console"""
    db_config_to_print = db_config.copy()
    if 'password' in db_config_to_print.keys():
        db_config_to_print['password'] = '*********'
    return db_config_to_print


def get_db_binding_config(db_config):
    """Removes name from db config because name is not needed for binding to pony orm database object"""
    db_binding_config = db_config.copy()
    if 'name' in db_binding_config.keys():
        db_binding_config.pop('name')
    return db_binding_config
