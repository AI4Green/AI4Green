import os
from typing import Literal, Union

import click
from compound_database import Compound_database_extraction as CDE
from compound_database import auxiliary, pubchem_download
from flask import current_app
from flask.cli import with_appcontext
from sources import models
from sources.extensions import db
from sqlalchemy import func


@click.command("drop-db")
@click.confirmation_option(prompt="Are you sure you want to drop the database?")
@with_appcontext
def drop_db() -> None:
    """
    Drop your database tables.
    """
    if not current_app.config["DEBUG"]:
        click.confirm(
            "You are about to drop the database in production mode. "
            "Are you sure you want to proceed?",
            abort=True,
        )
    db.reflect()
    db.drop_all()
    click.echo("Dropped the database.")


@click.command("seed-db")
@with_appcontext
def seed_db() -> None:
    """
    Seed your database with initial Pubchem data and configuration.
    """
    # unzip the download file, and populate the new database
    auxiliary.unzip_file(
        auxiliary.compound_database_dir / "Pubchem-Database",
        auxiliary.compound_database_dir / "Pubchem-Database_db.xml",
    )

    click.echo("Seeding element sustainability data")
    CDE.seed_element_sustainability_data()

    click.echo("Seeding hazard code data")
    CDE.seed_hazard_code_data()

    click.echo("Seeding user roles")
    CDE.seed_role_types()

    click.echo("Seeding compound data from pubchem download")
    compounds_file = os.path.join(
        auxiliary.compound_database_dir / "Pubchem-Database_db.xml"
    )
    CDE.extract_from_pubchem(compounds_file, current_app.config["COMPOUND_LIMIT"])

    click.echo("Seeding solvent data")
    CDE.seed_solvent_data()

    auxiliary.delete_temp_files()
    click.echo("Database seed completed.")
    db.session.commit()


@click.command("seed-users")
@with_appcontext
def seed_users() -> None:
    """
    Seeds predefined users, workgroups, and workbooks from `seed_data.yaml`
    """
    CDE.seed_predefined_data()
    click.echo("Seeded predefined data.")


@click.command("add-user")
@click.option("--username", prompt=True, type=str)
@click.option("--email", prompt=True, type=str)
@click.option("--fullname", prompt=True, type=str)
@click.password_option()
@click.option(
    "--role",
    prompt=True,
    type=click.Choice(["admin", "standard"], case_sensitive=False),
)
@with_appcontext
def add_user(
    username: str,
    email: str,
    fullname: str,
    password: str,
    role: Union[Literal["admin"], Literal["standard"]],
) -> None:
    """
    Adds a user to the database.

    Args:
        username: Username for the user.
        email: Email for the user.
        fullname: Full name of the user.
        password: Password for the user.
        role: Role to give the user.

    """
    role_model = (
        db.session.query(models.Role)
        .filter(func.lower(models.Role.name) == role)
        .first()
    )
    if role_model is None:
        raise ValueError(f"{role} not found, run `seed-db` first.")

    person = models.Person.create()
    models.User.create(
        username=username,
        email=email,
        fullname=fullname,
        person=person.id,
        password_hash=models.User.set_password(password),
        role=role_model.id,
    )
    click.echo(f"Created {username}")


@click.command("update-db")
@with_appcontext
def update_db() -> None:
    """
    Update your database with the latest Pubchem data.
    """
    list_of_updates = pubchem_download.find_database_updates()

    # Check for existing data
    check = db.session.query(models.Compound).first()
    if not check:
        click.echo("Cannot update empty database. Please create a new database.")
        return

    for update_date_link in list_of_updates:
        # Download the update.
        download_url = f"https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Monthly/{update_date_link}Extras/CID-LCSS.xml.gz"
        pubchem_download.download_pubchem_database(download_url, "Pubchem-update")

        # Unzip and make changes.
        auxiliary.unzip_file(
            auxiliary.compound_database_dir / "Pubchem-update",
            auxiliary.compound_database_dir / "Pubchem-update_db.xml",
        )
        compounds_file = os.path.join(
            auxiliary.compound_database_dir, "Pubchem-update_db.xml"
        )

        # Create tables in update
        db.create_all(bind_key="update")

        # Load data and update
        CDE.extract_from_pubchem(
            compounds_file, current_app.config["COMPOUND_LIMIT"], models.UpdateCompound
        )
        pubchem_download.apply_database_updates()

        db.drop_all(bind_key="update")
    auxiliary.delete_temp_files()
    auxiliary.write_date_to_file()
    click.echo("Database updated.")


@click.command("download-pubchem")
@with_appcontext
def download_pubchem() -> None:
    """
    This command will download the latest Pubchem database.
    """
    pubchem_download.database_download()
    click.echo("Pubchem download completed.")
