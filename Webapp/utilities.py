import os
import yaml
from shutil import copyfile

basedir = os.path.abspath(os.path.dirname(__file__))  # path to the user database
yaml_filepath = os.path.join(os.path.dirname(basedir), 'Webapp', 'configs.yaml')


def read_yaml(target_keys: list):
    """
    target_keys: iterates through keys to find nested data in yaml config
    """
    with open(yaml_filepath, "r") as stream:
        try:
            target_data = yaml.safe_load(stream)
            for key in target_keys:
                target_data = target_data[key]
        except yaml.YAMLError as exc:
            print(exc)
        return target_data


def read_yaml_db_config(db_name):
    db_config = read_yaml(["database_configurations", db_name])
    if db_config['provider'] == 'sqlite':
        if 'UNIT_TEST' in os.environ:
            copyfile(os.path.join(basedir, 'test_database.sqlite'), os.path.join(basedir, 'unit_test_database.sqlite'))
            db_config['filename'] = os.path.join(basedir, 'unit_test_database.sqlite')
        else:
            db_config['filename'] = os.environ.get('DATABASE_URL') or os.path.join(basedir, 'test_database.sqlite')
    return db_config



