"""
The code here downloads the PubChem database and the MD5 Checksum file.
It then confirms the successful download by running a checksum
"""
from auxiliary import *
from datetime import datetime
import Compound_database_extraction as CDE
import os
from pony.orm import select, commit, db_session
from utilities import read_yaml_db_config


def database_download():
    """Downloads the PubChem Database if it is not already present locally or if user requests a new version"""
    if os.path.isfile('Pubchem-Database'):
        time_downloaded = get_date('Download-Date.txt')
        print("\n\n")
        while True:
            # Ask the user whether to download new version or use existing present version
            redownload = input(
                f"A version of the PubChem database was downloaded on {time_downloaded}\nWould you like to "
                f"download the latest version?\nEnter 'y' to download the latest version or 'n' to use the "
                f"existing version.\n")
            if redownload in ['y', 'Y', 'n', 'N']:
                break
        if redownload in ['y', 'Y']:
            download_pubchem_database()
            write_date_to_file()
            return
        elif redownload == ['n', 'N']:
            print("Using existing database")
            return
    else:
        """
        If the database doesn't exist we download it
        """
        download_pubchem_database()
        write_date_to_file()
        return


def download_pubchem_database(download_url='http://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-LCSS.xml.gz',
                              save_file_as='Pubchem-Database'):
    """Download the database and verify the checksum"""
    print(f'Downloading PubChem database file: {download_url} and saving as {save_file_as}')
    download_file(download_url, save_file_as)
    download_file(download_url + '.md5', 'Checksum-new')
    checksum_file = os.path.join('Checksum-new')
    pub_chem_db_file = os.path.join(save_file_as)
    md5_integrity_check(pub_chem_db_file, checksum_file)


def find_database_updates():
    """
    To update the database we need to download the correct files
    and the checksums for them
    Returns a list of links from updates after the current database was downloaded
    """
    date_last_obj = get_date('Download-Date.txt')
    print(date_last_obj)
    list_of_links = find_files('https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Monthly/')
    refined_list_of_links = []
    for link in list_of_links:
        try:
            date_time_obj = datetime.strptime(link[:-1], '%Y-%m-%d')
            if date_last_obj <= date_time_obj:
                refined_list_of_links.append(link)
        except ValueError:
            pass
    print(f"Updates to apply: {refined_list_of_links}")
    return refined_list_of_links


def clear_update_table(db_update):
    """Clear the update table"""
    db_update.drop_all_tables(with_all_data=True)


def get_db_configs():
    """Loads db configurations from the yaml file."""
    configs_ls = [read_yaml_db_config("db_postgres"), read_yaml_db_config("db_sqlite")]
    return configs_ls


def databases_to_edit(db_config_ls_check):
    """User input for which databases to make/update"""
    db_config_ls = []
    # Determine whether to keep config in list or remove by asking for user input
    while True:
        make_config = input(f"Do you want to create a database using SQLite (recommended for new users) or "
                            f"Postgres (advanced)?\nEnter 's' for SQLite or 'p' for Postgres.\n")
        if make_config == "p":
            desired_provider = "postgres"
            break
        if make_config == "s":
            desired_provider = "sqlite"
            break

    for config in db_config_ls_check:
        if desired_provider == config["provider"]:
            db_config_ls.append(config)
            return db_config_ls


def create_or_update():
    """Ask the user whether to create or update a database"""
    while True:
        make_new_or_update = input(f"Would you like to create a new database or update the compound data in an "
                                   f"existing database?\nEnter 'c' to create a database "
                                   f"or 'u' to update the database.\n")
        if make_new_or_update in ['u', 'U', 'c', 'C']:
            break
    if make_new_or_update in ['c', 'C']:
        return 'create'
    elif make_new_or_update in ['u', 'U']:
        return 'update'


def create_database(db_config):
    """Handles database creation. Including old database destination, and populating the new database"""
    provider = db_config['provider']
    if 'filename' in db_config.keys():
        if os.path.isfile(db_config['filename']):
            determine_old_db_destination(db_config)
    # create database object
    db1 = CDE.open_database_compound(db_config)
    # unzip the download file, and populate the new database
    unzip_file('Pubchem-Database')
    Compounds_file = os.path.join('Pubchem-Database_db.xml')
    CDE.migrate_element_sustainability_data(db1)
    CDE.migrate_hazard_code_data(db1)
    CDE.add_role_types(db1)
    print("Extracting compound data from pubchem download")
    CDE.extract_from_pubchem(Compounds_file, db1, compound_limit)
    CDE.migrate_solvent_data(db1)
    print("All Compound data added")
    CDE.add_predefined_data(db1)


def update_database_handler(db_config):
    """Handles updating the database. Gets the list of updates and then applies each in turn"""
    print("This is where we update the existing database")
    list_of_updates = find_database_updates()
    # connect to the old database
    db_original = CDE.open_database_compound(db_config)
    with db_session:
        db_existing_data_check = select(x for x in db_original.Compound).first()
    if not db_existing_data_check:
        print("Cannot update empty database. Please create a new database.")
        exit()
    for update_date_link in list_of_updates:
        # make a temporary database for each update
        temp_update_db_config = {'provider': 'sqlite', 'filename': 'temp_update.sqlite', 'create_db': True}
        db_update = CDE.open_database_compound(temp_update_db_config)
        download_url = 'https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Monthly/' + update_date_link + 'Extras/CID-LCSS.xml.gz'
        download_pubchem_database(download_url, 'Pubchem-update')
        unzip_file('Pubchem-update')
        CDE.extract_from_pubchem('Pubchem-update_db.xml', db_update, compound_limit)
        apply_database_updates(db_original, db_update)
        clear_update_table(db_update)
        # write into date file down after completion
    delete_temp_files()
    print("Database updates complete")
    write_date_to_file()


@db_session
def apply_database_updates(db_original, db_update):
    """Compares the old database and the temporary database and applies the required updates"""
    compounds_in_old_db = set(select(old.CID for old in db_original.Compound)[:])
    compounds_in_update_db = set(select(new.CID for new in db_update.Compound)[:])
    compounds_in_both_dbs = compounds_in_update_db.intersection(compounds_in_old_db)
    update_db_entries(db_update, db_original, compounds_in_both_dbs)
    # delete old
    entries_to_be_deleted = compounds_in_old_db.difference(compounds_in_update_db)
    remove_old_entries(db_original, entries_to_be_deleted)
    # add new
    new_entries = compounds_in_update_db.difference(compounds_in_old_db)
    add_new_entries(db_original, db_update, new_entries)
    # check updates and then clear the temporary database
    check_db_updates(db_original, db_update, compounds_in_update_db)


@db_session
def update_db_entries(db_update, db_old, compounds_in_both_dbs):
    # if entries are the same dont update, if they differ, then update
    print(f"Compounds in both dbs to check for updates: {len(compounds_in_both_dbs)}")
    for cid in compounds_in_both_dbs:
        old_entry = select(old for old in db_old.Compound if old.CID == cid).first().to_dict()
        new_entry = select(new for new in db_update.Compound if new.CID == cid).first().to_dict()
        old_entry_id = old_entry.pop('id')
        del(new_entry['id'])
        if new_entry != old_entry:
            db_old.Compound[old_entry_id].set(**new_entry)
            commit()


@db_session
def remove_old_entries(db_old, entries_to_be_deleted):
    print(f"Removing {len(entries_to_be_deleted)} from database")
    for cid in entries_to_be_deleted:
        deletion_id = select(c for c in db_old.Compound if c.CID == cid).first().get_pk()
        db_old.Compound[deletion_id].delete()
        commit()


@db_session
def add_new_entries(db_old, db_update, entries_to_be_added):
    print(f"Adding {len(entries_to_be_added)} to database")
    for cid in entries_to_be_added:
        compound_to_add = select(new for new in db_update.Compound if new.CID == cid).first().to_dict()
        del (compound_to_add['id'])
        db_old.Compound(**compound_to_add)
        commit()


@db_session
def check_db_updates(db_old, db_update, compounds_in_db_update):
    for cid in compounds_in_db_update:
        # if i in cid_list_old:
        new_entry = select(new for new in db_update.Compound if new.CID == cid).first().to_dict()
        old_entry = select(old for old in db_old.Compound if old.CID == cid).first().to_dict()
        del(old_entry['id'], new_entry['id'])
        if new_entry != old_entry:
            print('Error with update to' + old_entry.get('CID'))


def determine_old_db_destination(db_config):
    """User input options for what to do with the existing database"""
    while True:
        existing_db_fate = input(
            f"File: {db_config['filename']} already exists.\nEnter 'x' to delete the existing database or 'a' to archive and make a new database."
            f"\n")
        if existing_db_fate in ['X', 'x', 'A', 'a']:
            break
    if existing_db_fate in ['X', 'x']:
        # delete db
        os.remove(db_config['filename'])
        print("Database deleted")
    if existing_db_fate in ['a', 'A']:
        # save to archive folder in compound database directory
        if not os.path.exists('archived_databases'):
            os.mkdir('archived_databases')
        archived_filename = os.path.abspath(
            os.path.join('archived_databases', 'db_' + datetime.today().strftime('%Y-%m-%d_%H-%M-%S') + '.sqlite'))
        os.rename(db_config['filename'], archived_filename)
        print("Database archived to: ", archived_filename)


if __name__ == '__main__':
    """
    To set up the postgres database locally follow these steps:
    1) Download the relevant version depending on OS: https://www.postgresql.org/download/
    2) Once installed navigate to the C:\Program Files\PostgreSQL\14\data\pg_hba.conf file
    3) Open the file with a text editor and change all 'METHOD' rows from 'scram-sha-256' to 'trust'
    source: https://stackoverflow.com/questions/64198359/pg-admin-4-password-for-postgres-user-when-trying-to-connect-to-postgresql-1
    4) reload the server config in postgres
    5) Run this python file to populate the database with the data required to run the AI4green web application.
    6) Use the same config in the main web app to connect to the now populated database to run AI4green
    
    To make a local database with a different path set the environmental variable "DB_FILEPATH".
    """
    # download database from pubchem or uses existing file
    database_download()
    # get database configurations
    db_config_ls_check = get_db_configs()
    # refine list through user input
    db_config_ls = databases_to_edit(db_config_ls_check)
    for db_config in db_config_ls:
        # for each db either create or update
        db_action = create_or_update()
        if db_action == 'create':
            create_database(db_config)
        elif db_action == 'update':
            update_database_handler(db_config)
