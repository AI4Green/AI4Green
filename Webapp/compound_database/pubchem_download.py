"""
The code here downloads the PubChem database and the MD5 Checksum file.
It then confirms the successful download by running a checksum
"""

import contextlib
import hashlib
import os
import sys
from datetime import datetime
from subprocess import run
from typing import Iterable, List, Literal, Optional, Union

import psycopg2
from compound_database.auxiliary import (
    compound_database_dir,
    download_file,
    find_files,
    get_date,
    md5_integrity_check,
    write_date_to_file,
)
from sources import models
from sources.extensions import db
from sqlalchemy.exc import IntegrityError

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))


def database_download():
    """Downloads the PubChem Database if it is not already present locally or if user requests a new version"""
    if os.path.isfile(compound_database_dir / "Pubchem-Database"):
        time_downloaded = get_date(compound_database_dir / "Download-Date.txt")
        print("\n\n")
        while True:
            # Ask the user whether to download new version or use existing present version
            redownload = input(
                f"A version of the PubChem database was downloaded on {time_downloaded}\nWould you like to "
                f"download the latest version?\nEnter 'y' to download the latest version or 'n' to use the "
                f"existing version.\n"
            )
            if redownload in ["y", "Y", "n", "N"]:
                break
        if redownload in ["y", "Y"]:
            download_pubchem_database()
            write_date_to_file()
            return
        elif redownload == ["n", "N"]:
            print("Using existing database")
            return
    else:
        """
        If the database doesn't exist we download it
        """
        download_pubchem_database()
        write_date_to_file()
        return


def download_pubchem_database(
    url: str = "http://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-LCSS.xml.gz",
    filename: str = "Pubchem-Database",
) -> None:
    """
    Download the database and verify the checksum.

    Args:
        url: The URL to download from.
        filename: Name of the file to save as.
    """
    print(f"Downloading PubChem database file: {url} and saving as {filename}")
    download_file(url, compound_database_dir / filename)
    download_file(f"{url}.md5", compound_database_dir / "Checksum-new")
    checksum_file = compound_database_dir / "Checksum-new"
    pub_chem_db_file = compound_database_dir / filename
    md5_integrity_check(pub_chem_db_file, checksum_file)


def find_database_updates() -> List[str]:
    """
    To update the database we need to download the correct files
    and the checksums for them.

    Returns:
        A list of links from updates after the current database was downloaded
    """
    date_last_obj = get_date(f"{compound_database_dir}/Download-Date.txt")
    print(f"Datebase last updated: {date_last_obj}")
    list_of_links = find_files("https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Monthly/")
    refined_list_of_links = []
    for link in list_of_links:
        with contextlib.suppress(ValueError):
            date_time_obj = datetime.strptime(link[:-1], "%Y-%m-%d")
            if date_last_obj <= date_time_obj:
                refined_list_of_links.append(link)
    print(f"Updates to apply: {refined_list_of_links}")
    return refined_list_of_links


def databases_to_edit(db_config_ls_check):
    """User input for which databases to make/update"""
    db_config_ls = []
    # Determine whether to keep config in list or remove by asking for user input
    while True:
        make_config = input(
            "Do you want to create a database using SQLite (recommended for new users) or "
            "Postgres (advanced)?\nEnter 's' for SQLite or 'p' for Postgres.\n"
        )
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


def create_or_update() -> Optional[Literal["create", "update"]]:
    """
    Ask the user whether to create or update a database.

    Returns:
        Create or Update based on the users CLI input.
    """
    while True:
        make_new_or_update = input(
            "Would you like to create a new database or update the compound data in an "
            "existing database?\nEnter 'c' to create a database "
            "or 'u' to update the database.\n"
        )
        if make_new_or_update in ["u", "U", "c", "C"]:
            break
    if make_new_or_update in ["c", "C"]:
        return "create"
    elif make_new_or_update in ["u", "U"]:
        return "update"


def determine_existing_db_destination(db_config):
    """Either sqlite or postgres method for handling existing database"""
    if db_config["provider"] == "sqlite":
        determine_existing_db_destination_sqlite(db_config)
    elif db_config["provider"] == "postgres":
        determine_existing_db_destination_postgres(db_config)
    else:
        exit("Only 'sqlite' and 'postgres' database providers are supported")


def determine_existing_db_destination_sqlite(db_config):
    """Archives or removes existing sqlite database if present"""
    if check_if_existing_sqlite_db(db_config) is True:
        existing_db_fate = determine_db_fate()
        if existing_db_fate in ["X", "x"]:
            # delete db
            os.remove(db_config["filename"])
            print("Database deleted")
        if existing_db_fate in ["a", "A"]:
            # save to archive folder in compound database directory
            archived_filename = get_archive_file(".sqlite")
            os.rename(db_config["filename"], archived_filename)
            print("Database archived to: ", archived_filename)


def determine_existing_db_destination_postgres(db_config):
    """Archives or removes existing postgres database is present"""
    cursor = connect_to_postgres_db(db_config)
    if check_if_existing_postgres_db(cursor, db_config):
        existing_db_fate = determine_db_fate()
        if existing_db_fate in ["x", "X"]:
            # delete db
            postgres_delete_all_data(cursor)
        if existing_db_fate in ["a", "A"]:
            # use pg dump on the cli to save the sql file to the archive folder in compound database directory
            archived_filename = get_archive_file(".sql")
            pg_dump_to_sql(db_config, archived_filename)
            # cascade delete all data and then check tables have been dropped
            postgres_delete_all_data(cursor)


def check_if_existing_postgres_db(cursor, db_config):
    """Returns True if existing postgres database is found"""
    # checks if database exists in the list of databases
    cursor.execute("SELECT datname FROM pg_database;")
    list_of_databases = cursor.fetchall()
    # check database exists in postgres
    if (f'{db_config["database"]}',) in list_of_databases:
        cursor.execute(
            """SELECT table_name FROM information_schema.tables WHERE table_schema = 'public'"""
        )
        database_tables = cursor.fetchall()
        # check if any data is in the database - 'Element' is the first table added to the database
        if ("Element",) in database_tables:
            print("Database present")
            return True
    print("Database not present")
    return False


def check_if_existing_sqlite_db(db_config):
    """Returns True is existing sqlite database is found"""
    if os.path.isfile(compound_database_dir / db_config["filename"]):
        return True
    return False


def get_archive_file(file_extension):
    """Returns a filepath to archive the database file and makes the directory if needed"""
    if not os.path.exists(compound_database_dir / "archived_databases"):
        os.mkdir(compound_database_dir / "archived_databases")
    file_path = (
        compound_database_dir
        / "archived_databases"
        / datetime.today().strftime("%Y-%m-%d_%H-%M-%S")
    )
    return file_path.with_suffix(file_extension)


def determine_db_fate():
    """Returns the user input for decision on whether to delete or archive existing database"""
    while True:
        existing_db_fate = input(
            "\nDatabase with this configuration already exists.\nEnter 'x' to delete the existing database or 'a' to archive and make a new database."
            "\n"
        )
        if existing_db_fate in ["X", "x", "A", "a"]:
            break
    return existing_db_fate


def connect_to_postgres_db(db_config):
    """Connects to the postgres database and returns a cursor object"""
    connection = psycopg2.connect(
        database=db_config["database"],
        user=db_config["user"],
        password=db_config["password"],
        host=db_config["host"],
        port=db_config["port"],
    )
    connection.autocommit = True
    return connection.cursor()


def postgres_delete_all_data(cursor):
    """Does a cascade delete of the current schema to remove all tables. Then makes a new blank schema ready to accept
    any changes to the data format."""
    cursor.execute("DROP SCHEMA public CASCADE")
    cursor.execute("CREATE SCHEMA public;")
    # checks that there no tables left in the postgres database
    assert (
        cursor.execute("select table_schema, table_name from information_schema.tables")
        is None
    )


def pg_dump_to_sql(db_config, archived_filename):
    """Uses a terminal command to run pg_dump to make an archive of an old postgres database"""
    command = (
        f'pg_dump --host={db_config["host"]} --username={db_config["user"]} '
        f"--dbname=postgres --file={archived_filename}"
    )
    command_result = run(command, shell=True, env={"PGPASSWORD": db_config["password"]})
    # check pg dump was successful before deleting
    assert (
        command_result.returncode == 0
        and os.path.exists(compound_database_dir / archived_filename) is True
    )


def apply_database_updates() -> None:
    """
    Compares the main database and the temporary database and applies the required updates.
    """
    compounds_in_old_db = {
        row for (row,) in db.session.query(models.Compound.cid).all()
    }
    compounds_in_update_db = {
        row for (row,) in db.session.query(models.UpdateCompound.cid).all()
    }

    # 1. Delete old
    ids_to_delete = compounds_in_old_db.difference(compounds_in_update_db)
    compounds_to_delete = (
        db.session.query(models.Compound)
        .filter(models.Compound.cid.in_(ids_to_delete))
        .all()
    )
    remove_old_entries(compounds_to_delete)

    # 2. Update existing
    ids_in_both_dbs = compounds_in_old_db.intersection(compounds_in_update_db)
    compounds_to_update = (
        db.session.query(models.Compound)
        .filter(models.Compound.cid.in_(ids_in_both_dbs))
        .all()
    )
    update_db_entries(compounds_to_update)

    # 3. Add new
    ids_to_add = compounds_in_update_db.difference(compounds_in_old_db)
    compounds_to_add = (
        db.session.query(models.UpdateCompound)
        .filter(models.UpdateCompound.cid.in_(ids_to_add))
        .all()
    )
    add_new_entries(compounds_to_add)

    # 4. Check updates
    check_db_updates(compounds_in_update_db)


def update_db_entries(compounds: Iterable[models.Compound]) -> None:
    """
    Updates a list of compounds.
    If a compound is the same do not update, if they differ, then update.

    Args:
        compounds: Compounds to update that exist in both databases.
    """
    print(f"Compounds in both dbs to check for updates: {len(compounds)}")
    excluded_fields = ["id", "solvent", "error_report"]

    for c in compounds:
        new = (
            db.session.query(models.UpdateCompound)
            .filter(models.UpdateCompound.cid == c.cid)
            .first()
        )

        c_hash = compute_compound_hash(c)
        new_hash = compute_compound_hash(new)
        if c_hash != new_hash:
            for field in c.__table__.columns:
                if field.name not in excluded_fields:
                    field_name = field.name
                    new_value = getattr(new, field_name)
                    setattr(c, field_name, new_value)
    db.session.commit()


def remove_old_entries(compounds: Iterable[models.Compound]) -> None:
    """
    Deletes a list of compounds in the database.
    TODO: Change this to bulk delete.

    Args:
        compounds: the iterable of compounds to delete.
    """
    print(f"Removing {len(compounds)} from database")
    for compound in compounds:
        db.session.delete(compound)
    db.session.commit()


def add_new_entries(compounds: Iterable[models.UpdateCompound]) -> None:
    """
    Adds new compounds to the database.

    Args:
        compounds: the list of compounds to add.
    """
    print(f"Adding {len(compounds)} to database")
    for compound in compounds:
        compound_data = compound.to_dict()
        new_compound = models.Compound(**compound_data)
        db.session.add(new_compound)
        try:
            db.session.commit()
        except IntegrityError:
            db.session.rollback()
            print(f"Duplicate CAS key found for compound CAS: {new_compound.cas}.")


def check_db_updates(compounds: Iterable[str]) -> None:
    """
    Checks two databases have updated, using the list of compound ids to check against.

    Args:
        compounds: the list of Compounds ids to check.
    """
    for c in compounds:
        new_entry = (
            db.session.query(models.UpdateCompound)
            .filter(models.UpdateCompound.cid == c)
            .first()
        )
        old_entry = (
            db.session.query(models.Compound).filter(models.Compound.cid == c).first()
        )

        if new_entry is None or old_entry is None:
            print(f"Entry missing: {c}")
            continue

        if compute_compound_hash(new_entry) != compute_compound_hash(old_entry):
            print(f"Error with update to {old_entry.cid}")


def compute_compound_hash(
    compound: Union[models.Compound, models.UpdateCompound]
) -> str:
    """
    Computes a hash value for a compound based on its field values.

    Args:
        compound: The compound object.

    Returns:
        A string representing the hash value.
    """
    excluded_fields = ["id", "solvent", "error_report"]

    # Combine the field values into a single string and compute the hash
    compound_data = "".join(
        str(getattr(compound, field.name))
        for field in compound.__table__.columns
        if field.name not in excluded_fields
    )
    return hashlib.md5(compound_data.encode()).hexdigest()
