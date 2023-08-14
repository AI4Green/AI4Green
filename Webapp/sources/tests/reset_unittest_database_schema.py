from utilities import read_yaml_db_config, read_yaml


def reset_unittest_db(db_config):
    """
    Deletes and then remakes the database
    """
    from compound_database.pubchem_download import create_database,\
        connect_to_postgres_db, postgres_delete_all_data
    cursor = connect_to_postgres_db(db_config)
    while True:
        delete_unitttest_db = input('\nUnit test database schema does not match models. Enter "X" to delete and remake '
                                    'an updated version of the unit test database or "E" to exit this script\n')
        if delete_unitttest_db in ['E', 'e']:
            exit()
        if delete_unitttest_db in ['X', 'x']:
            postgres_delete_all_data(cursor)
            create_database(db_config)
            break


if __name__ == '__main__':
    """
    Run this file to reset the unittest database
    """
    database_type = read_yaml(["database_configurations", "active_db"])
    db_config = read_yaml_db_config(database_type)
    db_config['database'] = 'unit_test'
    reset_unittest_db(db_config)