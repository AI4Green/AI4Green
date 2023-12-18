from datetime import datetime
from typing import Optional, Tuple

from flask import abort
from flask_login import current_user
from sources import models, services
from sources.extensions import db
from sqlalchemy import func


def get_smiles(primary_key: Tuple[str, int]) -> str:
    """
    Gets the novel compound's SMILES string from the primary key if the entry has the SMILES attribute

    Returns:
        The SMILES string corresponding to the primary key or None
    """
    primary_key = (primary_key[0], int(primary_key[1]))
    workbook = services.workbook.get(primary_key[1])
    if current_user.Person not in workbook.users:
        abort(401)

    return (
        db.session.query(models.NovelCompound.smiles)
        .filter(models.NovelCompound.name == primary_key[0])
        .join(models.WorkBook)
        .filter(models.WorkBook.id == primary_key[1])
        .first()
    )[0]


def get_novel_compound_from_name_and_workbook(
    name: str, workbook: models.WorkBook
) -> models.NovelCompound:
    """
    Retrieves a novel compound by name and workbook.

    Args:
        name: Compound name.
        workbook: Workbook model.

    Returns:
        NovelCompound model.
    """
    return (
        db.session.query(models.NovelCompound)
        .filter(func.lower(models.NovelCompound.name) == name.lower())
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook.id)
        .first()
    )


def get_novel_compound_from_inchi_and_workbook(
    inchi: str, workbook: models.WorkBook
) -> models.NovelCompound:
    """
    Retrieves a novel compound by name and workbook.

    Args:
        inchi: Compound inchi.
        workbook: Workbook model.

    Returns:
        NovelCompound model.
    """
    return (
        db.session.query(models.NovelCompound)
        .filter(models.NovelCompound.inchi == inchi)
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook.id)
        .first()
    )


def get_novel_compound_from_cas_and_workbook(
    cas: str, workbook: models.WorkBook
) -> models.NovelCompound:
    """
    Retrieves a novel compound by CAS and workbook.

    Args:
        cas: CAS number.
        workbook: Workbook model.

    Returns:
        NovelCompound model.
    """
    return (
        db.session.query(models.NovelCompound)
        .filter(models.NovelCompound.cas == cas)
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook.id)
        .first()
    )


def add(
    name: str,
    cas: Optional[str],
    mol_formula: str,
    mol_weight: float,
    density: float,
    concentration: float,
    hazards: str,
    smiles: str,
    inchi: Optional[str],
    inchi_key: Optional[str],
    workbook_id: int,
) -> models.NovelCompound:
    """
    Creates a novel compound in the database.

    Args:
        name: Compound name.
        cas: CAS number.
        mol_formula: Molecule formula.
        mol_weight: Molecular weight.
        density: Density.
        concentration: Concentration.
        hazards: Hazard codes.
        smiles: SMILES notation.
        inchi: InChI notation.
        inchi_key: InChIKey.
        workbook_id: Workbook ID.
        current_time: Current timestamp.

    Returns:
        NovelCompound model.
    """
    nc = models.NovelCompound(
        name=name,
        cas=cas,
        molec_formula=mol_formula,
        molec_weight=mol_weight,
        density=density,
        concentration=concentration,
        hphrase=hazards,
        smiles=smiles,
        inchi=inchi,
        inchikey=inchi_key,
        workbook=workbook_id,
    )
    db.session.add(nc)
    db.session.commit()
    return nc
