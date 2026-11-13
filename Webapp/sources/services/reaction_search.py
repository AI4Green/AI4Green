"""
Searching the reactions a user can see, for the React search page.

A search combines any of:
- a keyword, matched against the reaction's name, description and ID
- a reagent or solvent name
- a creation date range
- a drawn structure, matched exactly or as a substructure against the
  reaction's reactants, reagents, solvents and/or products

Security notes:
- Only active reactions in workbooks the current user belongs to are searched.
  A workgroup or workbook the user can't see raises ReactionSearchScopeError,
  which the API turns into a 404 so it doesn't reveal what exists.
- Input is validated here. Problems raise ReactionSearchValidationError (or
  StructureSearchError) with a message that is safe to show to the user.
"""

import json
from dataclasses import dataclass, field
from datetime import date, datetime, time, timedelta
from typing import Any, Dict, FrozenSet, Iterable, List, Optional, Set

from sources import models, services
from sources.auxiliary import get_workbooks, get_workgroups
from sources.extensions import db
from sources.models.reaction import ReactionType
from sources.services.structure_search import StructureQuery
from sqlalchemy import func, or_

# The parts of a reaction a structure can be matched against
ROLES = ("reactant", "reagent", "solvent", "product")
MAX_RESULTS = 200
MAX_TEXT_FILTER_LENGTH = 200
DESCRIPTION_PREVIEW_LENGTH = 250
# Longest novel compound primary key string, matching
# services.novel_compound.reform_novel_compound_primary_key
MAX_NOVEL_PRIMARY_KEY_LENGTH = 350


class ReactionSearchValidationError(ValueError):
    """The search request is invalid. The message is safe to show to the user."""


class ReactionSearchScopeError(LookupError):
    """The requested workgroup or workbook doesn't exist or the user can't see it."""


@dataclass(frozen=True)
class ReactionSearchCriteria:
    """
    What to search for. Empty fields don't filter.

    Attributes:
        workgroup: Only search this workgroup, or every workgroup if None.
        workbook: Only search this workbook of `workgroup`, or all of them if None.
        keyword: Text to find in the name, description or reaction ID.
        reagent: Text to find in a reagent or solvent name.
        date_from: Earliest creation date, inclusive.
        date_to: Latest creation date, inclusive.
        structure: The drawn structure to match.
        roles: Which parts of the reaction the structure is matched against.
    """

    workgroup: Optional[str] = None
    workbook: Optional[str] = None
    keyword: Optional[str] = None
    reagent: Optional[str] = None
    date_from: Optional[date] = None
    date_to: Optional[date] = None
    structure: Optional[StructureQuery] = None
    roles: FrozenSet[str] = frozenset(ROLES)


@dataclass
class ReactionSearchResult:
    """A matching reaction and the roles the structure was found in (if any)."""

    reaction: models.Reaction
    matched_roles: List[str] = field(default_factory=list)


@dataclass
class ReactionSearchResults:
    """The matches, and whether there were more than MAX_RESULTS."""

    results: List[ReactionSearchResult]
    truncated: bool


# --- Reading the request ----------------------------------------------------


def _read_text(value: Any, label: str) -> Optional[str]:
    """
    Reads an optional text filter.

    Args:
        value: The value from the request.
        label: The field's name for error messages.

    Returns:
        The stripped text, or None if it's empty.
    """
    if value is None:
        return None
    if not isinstance(value, str):
        raise ReactionSearchValidationError(f"{label} must be text.")
    value = value.strip()
    if len(value) > MAX_TEXT_FILTER_LENGTH:
        raise ReactionSearchValidationError(
            f"{label} must be at most {MAX_TEXT_FILTER_LENGTH} characters."
        )
    return value or None


def _read_date(value: Any, label: str) -> Optional[date]:
    """
    Reads an optional YYYY-MM-DD date.

    Args:
        value: The value from the request.
        label: The field's name for error messages.

    Returns:
        The date, or None if it's empty.
    """
    if value in (None, ""):
        return None
    if not isinstance(value, str):
        raise ReactionSearchValidationError(f"{label} must be a date.")
    try:
        return date.fromisoformat(value)
    except ValueError:
        raise ReactionSearchValidationError(f"{label} must be a date.") from None


def _read_roles(value: Any) -> FrozenSet[str]:
    """
    Reads which parts of the reaction to match the structure against.

    Args:
        value: A list of role names, or None for all of them.

    Returns:
        The roles.
    """
    if value is None:
        return frozenset(ROLES)
    if not isinstance(value, list) or not all(role in ROLES for role in value):
        raise ReactionSearchValidationError(
            "Choose where to match the structure: reactants, reagents, solvents or products."
        )
    if not value:
        raise ReactionSearchValidationError(
            "Choose at least one place to match the structure."
        )
    return frozenset(value)


def read_criteria(data: Dict[str, Any]) -> ReactionSearchCriteria:
    """
    Builds search criteria from the JSON request body.

    Args:
        data: The body: `workgroup`, `workbook`, `keyword`, `reagent`, `dateFrom`,
            `dateTo` and `structure` (`{smiles, mode, roles}`), all optional.

    Returns:
        The criteria.

    Raises:
        ReactionSearchValidationError: If a field is invalid or nothing was searched for.
        StructureSearchError: If the drawn structure can't be searched.
    """
    workgroup = _read_text(data.get("workgroup"), "Workgroup")
    workbook = _read_text(data.get("workbook"), "Workbook")
    if workbook and not workgroup:
        raise ReactionSearchValidationError("Choose a workgroup to pick a workbook.")

    date_from = _read_date(data.get("dateFrom"), "The start date")
    date_to = _read_date(data.get("dateTo"), "The end date")
    if date_from and date_to and date_from > date_to:
        raise ReactionSearchValidationError(
            "The start date must be on or before the end date."
        )

    structure = None
    roles = frozenset(ROLES)
    structure_data = data.get("structure")
    if structure_data is not None:
        if not isinstance(structure_data, dict):
            raise ReactionSearchValidationError("The structure is invalid.")
        mode = services.structure_search.read_mode(structure_data.get("mode"))
        structure = services.structure_search.build_query(
            structure_data.get("smiles"), mode
        )
        roles = _read_roles(structure_data.get("roles"))

    criteria = ReactionSearchCriteria(
        workgroup=workgroup,
        workbook=workbook,
        keyword=_read_text(data.get("keyword"), "The keyword"),
        reagent=_read_text(data.get("reagent"), "The reagent name"),
        date_from=date_from,
        date_to=date_to,
        structure=structure,
        roles=roles,
    )
    if not any(
        (
            criteria.keyword,
            criteria.reagent,
            criteria.date_from,
            criteria.date_to,
            criteria.structure,
        )
    ):
        raise ReactionSearchValidationError(
            "Enter a keyword, reagent, date or structure to search for."
        )
    return criteria


# --- Scope -------------------------------------------------------------------


def accessible_workbooks() -> List[Dict[str, Any]]:
    """
    The workgroups the current user belongs to, each with the workbooks they
    belong to, for the search filters.

    Returns:
        `[{name, workbooks: [name, ...]}, ...]`, sorted by workgroup name.
    """
    return [
        {"name": workgroup, "workbooks": sorted(get_workbooks(workgroup))}
        for workgroup in get_workgroups()
    ]


def _workbook_ids_in_scope(
    workgroup: Optional[str], workbook: Optional[str]
) -> List[int]:
    """
    The ids of the workbooks to search, limited to those the user belongs to.

    Workbook names are only unique within a workgroup, so workbooks are matched
    by (workgroup, workbook) pair rather than by name alone.

    Args:
        workgroup: The chosen workgroup, or None for all of them.
        workbook: The chosen workbook, or None for all of them.

    Returns:
        The workbook ids.

    Raises:
        ReactionSearchScopeError: If the chosen workgroup or workbook isn't one
            the user belongs to.
    """
    user_workgroups = get_workgroups()
    if workgroup is not None and workgroup not in user_workgroups:
        raise ReactionSearchScopeError()
    workgroups = [workgroup] if workgroup is not None else user_workgroups

    allowed_pairs: Set[tuple] = {
        (workgroup_name, workbook_name)
        for workgroup_name in workgroups
        for workbook_name in get_workbooks(workgroup_name)
        if workbook is None or workbook_name == workbook
    }
    if workbook is not None and not allowed_pairs:
        raise ReactionSearchScopeError()
    if not allowed_pairs:
        return []

    rows = (
        db.session.query(
            models.WorkBook.id, models.WorkBook.name, models.WorkGroup.name
        )
        .join(models.WorkGroup, models.WorkBook.group == models.WorkGroup.id)
        .filter(models.WorkGroup.name.in_(workgroups))
        .all()
    )
    return [
        workbook_id
        for workbook_id, workbook_name, workgroup_name in rows
        if (workgroup_name, workbook_name) in allowed_pairs
    ]


# --- Filtering ---------------------------------------------------------------


def _like_pattern(text: str) -> str:
    """
    Makes a case-insensitive "contains" pattern, escaping LIKE wildcards so
    the user's text is matched literally.

    Args:
        text: The user's text.

    Returns:
        The pattern, for use with `escape="\\"`.
    """
    escaped = text.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_")
    return f"%{escaped}%"


def _candidate_reactions(
    criteria: ReactionSearchCriteria, workbook_ids: List[int]
) -> List[models.Reaction]:
    """
    The active reactions in scope that pass the filters the database can apply:
    keyword and creation dates. Sorted by name.

    Args:
        criteria: The search criteria.
        workbook_ids: The workbooks to search.

    Returns:
        The candidate reactions.
    """
    query = db.session.query(models.Reaction).filter(
        models.Reaction.status == "active",
        models.Reaction.workbooks.in_(workbook_ids),
    )
    if criteria.keyword:
        pattern = _like_pattern(criteria.keyword)
        query = query.filter(
            or_(
                models.Reaction.name.ilike(pattern, escape="\\"),
                models.Reaction.description.ilike(pattern, escape="\\"),
                models.Reaction.reaction_id.ilike(pattern, escape="\\"),
            )
        )
    if criteria.date_from:
        query = query.filter(
            models.Reaction.time_of_creation
            >= datetime.combine(criteria.date_from, time.min)
        )
    if criteria.date_to:
        # inclusive: anything before the start of the next day
        query = query.filter(
            models.Reaction.time_of_creation
            < datetime.combine(criteria.date_to + timedelta(days=1), time.min)
        )
    return query.order_by(func.lower(models.Reaction.name), models.Reaction.id).all()


def _reaction_table_data(reaction: models.Reaction) -> Dict[str, Any]:
    """
    The reaction table data as a dict. It's saved with json.dumps into a JSON
    column, so it usually comes back as a string.

    Args:
        reaction: The reaction.

    Returns:
        The data, or an empty dict if it can't be read.
    """
    data = reaction.reaction_table_data
    if isinstance(data, str):
        try:
            data = json.loads(data)
        except ValueError:
            return {}
    return data if isinstance(data, dict) else {}


def _has_reagent_named(reaction: models.Reaction, text: str) -> bool:
    """
    Whether any reagent or solvent name contains the text, ignoring case.

    Args:
        reaction: The reaction.
        text: The name to look for.

    Returns:
        True if a reagent or solvent name contains it.
    """
    data = _reaction_table_data(reaction)
    names = list(data.get("reagent_names") or []) + list(
        data.get("solvent_names") or []
    )
    needle = text.casefold()
    return any(isinstance(name, str) and needle in name.casefold() for name in names)


def _solvent_smiles_by_key(
    reactions: Iterable[models.Reaction], workbook_ids: List[int]
) -> Dict[str, str]:
    """
    Looks up the SMILES of every solvent used by the reactions, in as few
    queries as possible.

    Reaction.solvent holds primary keys as strings: an integer for a Compound,
    or "('name', workbook_id)" for a NovelCompound. services.all_compounds.
    get_smiles_list resolves these, but it makes one query per key, raises on
    keys that no longer exist, and aborts the request for novel compounds in
    other workbooks, so it isn't suited to scanning many reactions. Instead,
    this batches the compound lookups, reuses the novel compound key parser,
    and only reads novel compounds from workbooks being searched.

    Args:
        reactions: The reactions whose solvents are needed.
        workbook_ids: The workbooks the user may see.

    Returns:
        SMILES by primary key string. Keys that can't be resolved are left out.
    """
    keys = {
        str(key) for reaction in reactions for key in (reaction.solvent or []) if key
    }
    smiles_by_key: Dict[str, str] = {}

    compound_ids = {int(key) for key in keys if key.isdigit()}
    if compound_ids:
        rows = (
            db.session.query(models.Compound.id, models.Compound.smiles)
            .filter(models.Compound.id.in_(compound_ids))
            .all()
        )
        smiles_by_key.update({str(id_): smiles for id_, smiles in rows if smiles})

    allowed_workbooks = set(workbook_ids)
    for key in keys:
        if key.isdigit() or len(key) > MAX_NOVEL_PRIMARY_KEY_LENGTH:
            continue
        try:
            (
                name,
                workbook_id,
            ) = services.novel_compound.reform_novel_compound_primary_key(key)
        except (AttributeError, ValueError):
            # not a novel compound key, e.g. a placeholder value
            continue
        if not name or workbook_id not in allowed_workbooks:
            continue
        smiles = (
            db.session.query(models.NovelCompound.smiles)
            .filter(
                func.lower(models.NovelCompound.name) == name.lower(),
                models.NovelCompound.workbook == workbook_id,
            )
            .scalar()
        )
        if smiles:
            smiles_by_key[key] = smiles
    return smiles_by_key


def _split_reaction_smiles(reaction_smiles: Optional[str]) -> tuple:
    """
    Splits reaction SMILES into reactant and product SMILES lists.

    Args:
        reaction_smiles: e.g. "C.CC>>CP".

    Returns:
        (reactants, products), each empty if the SMILES isn't a reaction.
    """
    if not reaction_smiles or ">>" not in reaction_smiles:
        return [], []
    reactants, products = reaction_smiles.split(">>", 1)
    return reactants.split("."), products.split(".")


def _flatten_polymer_components(smiles_list: Iterable[str]) -> List[str]:
    """
    Polymer reactions store some components as JSON lists of SMILES. This
    flattens them into one list, leaving plain SMILES as they are.

    Args:
        smiles_list: Stored SMILES, some of which may be JSON lists.

    Returns:
        A flat list of SMILES.
    """
    flat: List[str] = []
    for smiles in smiles_list:
        try:
            loaded = json.loads(smiles)
        except (TypeError, ValueError):
            flat.append(smiles)
            continue
        if isinstance(loaded, list):
            flat.extend(item for item in loaded if isinstance(item, str))
        else:
            flat.append(smiles)
    return flat


def _components_by_role(
    reaction: models.Reaction, solvent_smiles_by_key: Dict[str, str]
) -> Dict[str, List[str]]:
    """
    Every SMILES in the reaction, grouped by role.

    Reactants and products come from both the reaction SMILES and the saved
    lists, since either can be missing something the other has.

    Args:
        reaction: The reaction.
        solvent_smiles_by_key: Solvent SMILES from _solvent_smiles_by_key.

    Returns:
        SMILES lists keyed by role.
    """
    reactants_from_smiles, products_from_smiles = _split_reaction_smiles(
        reaction.reaction_smiles
    )
    reactants = list(reaction.reactants or [])
    products = list(reaction.products or [])
    if reaction.reaction_type == ReactionType.POLYMER:
        reactants = _flatten_polymer_components(reactants)
        products = _flatten_polymer_components(products)

    return {
        "reactant": list(set(reactants_from_smiles + reactants)),
        "reagent": list(reaction.reagents or []),
        "solvent": [
            solvent_smiles_by_key[str(key)]
            for key in (reaction.solvent or [])
            if str(key) in solvent_smiles_by_key
        ],
        "product": list(set(products_from_smiles + products)),
    }


def _matched_roles(
    reaction: models.Reaction,
    structure: StructureQuery,
    roles: FrozenSet[str],
    solvent_smiles_by_key: Dict[str, str],
) -> List[str]:
    """
    The roles in which the reaction contains the searched structure.

    Args:
        reaction: The reaction.
        structure: The prepared structure query.
        roles: The roles to look in.
        solvent_smiles_by_key: Solvent SMILES from _solvent_smiles_by_key.

    Returns:
        The matching roles, in ROLES order. Empty if there's no match.
    """
    # polymer repeat units are only stored on polymer reactions
    if structure.is_polymer and reaction.reaction_type != ReactionType.POLYMER:
        return []
    components = _components_by_role(reaction, solvent_smiles_by_key)
    return [
        role
        for role in ROLES
        if role in roles
        and services.structure_search.any_match(structure, components[role])
    ]


# --- Searching ---------------------------------------------------------------


def search_reactions(criteria: ReactionSearchCriteria) -> ReactionSearchResults:
    """
    Finds the reactions matching every given criterion, sorted by name.

    Args:
        criteria: What to search for.

    Returns:
        Up to MAX_RESULTS matches, and whether more were found.

    Raises:
        ReactionSearchScopeError: If the chosen workgroup or workbook isn't one
            the user belongs to.
    """
    workbook_ids = _workbook_ids_in_scope(criteria.workgroup, criteria.workbook)
    if not workbook_ids:
        return ReactionSearchResults(results=[], truncated=False)

    candidates = _candidate_reactions(criteria, workbook_ids)
    if criteria.reagent:
        candidates = [
            reaction
            for reaction in candidates
            if _has_reagent_named(reaction, criteria.reagent)
        ]

    if criteria.structure is None:
        return ReactionSearchResults(
            results=[
                ReactionSearchResult(reaction) for reaction in candidates[:MAX_RESULTS]
            ],
            truncated=len(candidates) > MAX_RESULTS,
        )

    solvent_smiles_by_key = (
        _solvent_smiles_by_key(candidates, workbook_ids)
        if "solvent" in criteria.roles
        else {}
    )
    results: List[ReactionSearchResult] = []
    for reaction in candidates:
        matched_roles = _matched_roles(
            reaction, criteria.structure, criteria.roles, solvent_smiles_by_key
        )
        if not matched_roles:
            continue
        if len(results) == MAX_RESULTS:
            return ReactionSearchResults(results=results, truncated=True)
        results.append(ReactionSearchResult(reaction, matched_roles))
    return ReactionSearchResults(results=results, truncated=False)


# --- Serialising -------------------------------------------------------------


def _iso(value: Optional[datetime]) -> Optional[str]:
    """
    Formats a datetime for JSON.

    Args:
        value: The datetime, or None.

    Returns:
        The ISO 8601 string, or None.
    """
    return value.isoformat() if value else None


def _safe_image(image: Optional[str]) -> Optional[str]:
    """
    The reaction image, only if it's an image data URL. The front end puts it
    in an <img src>, so anything else (e.g. a javascript: URL) is dropped.

    Args:
        image: The stored reaction image.

    Returns:
        The data URL, or None.
    """
    if isinstance(image, str) and image.startswith("data:image/"):
        return image
    return None


def serialise_result(result: ReactionSearchResult) -> Dict[str, Any]:
    """
    A search result for the front end.

    services.reaction.to_dict isn't reused: it builds the dicts for the Jinja
    reaction list, so it also loads each reaction's addenda and full table data
    and formats dates with str() rather than ISO 8601.

    Args:
        result: The matching reaction and the roles it matched in.

    Returns:
        The reaction's summary, in camelCase.
    """
    reaction = result.reaction
    creator_user = reaction.creator_person.user if reaction.creator_person else None
    description = reaction.description or ""
    if len(description) > DESCRIPTION_PREVIEW_LENGTH:
        description = description[: DESCRIPTION_PREVIEW_LENGTH - 1] + "..."
    return {
        "id": reaction.id,
        "reactionId": reaction.reaction_id,
        "name": reaction.name,
        "description": description,
        "workgroup": reaction.workbook.WorkGroup.name,
        "workbook": reaction.workbook.name,
        "creatorUsername": creator_user.username if creator_user else None,
        "timeOfCreation": _iso(reaction.time_of_creation),
        "timeOfUpdate": _iso(reaction.time_of_update),
        "completionStatus": reaction.complete,
        "image": _safe_image(reaction.reaction_image),
        "matchedRoles": result.matched_roles,
    }


def serialise_results(results: ReactionSearchResults) -> Dict[str, Any]:
    """
    The search results for the front end.

    Args:
        results: The search results.

    Returns:
        `{results: [...], truncated}`.
    """
    return {
        "results": [serialise_result(result) for result in results.results],
        "truncated": results.truncated,
    }
