"""
Journeys: creating a reaction in a workbook, building its reaction table from
a sketch, and deleting it.

Reactions are drawn with compounds seeded by tests/populate_database.py, so the
reaction table can name them without looking anything up online.
"""

import re

from e2e.pages.reaction_constructor_page import ReactionConstructorPage
from e2e.pages.reaction_list_page import ReactionListPage
from e2e.seed_data import SEEDED_REACTION_NAME, TEST_WORKBOOK, TEST_WORKGROUP
from playwright.sync_api import Page, expect

# methane + ethane >> methylphosphine, seeded as these compounds
REACTION_SMILES = "C.CC>>CP"
REACTANT_NAMES = ["Testoic Acid", "Testamine"]
PRODUCT_NAME = "Testproductine"


def test_user_creates_reaction(logged_in_page: Page) -> None:
    """A new reaction opens in the sketcher, and is listed in its workbook."""
    reaction_name = "E2E Created Reaction"

    reaction_list = ReactionListPage(logged_in_page)
    reaction_list.goto(TEST_WORKGROUP, TEST_WORKBOOK)
    reaction_list.create_reaction(reaction_name)

    reaction_page = ReactionConstructorPage(logged_in_page)
    expect(logged_in_page).to_have_url(re.compile(r"/sketcher/"))
    expect(reaction_page.heading).to_contain_text(reaction_name)

    reaction_list.goto(TEST_WORKGROUP, TEST_WORKBOOK)
    expect(reaction_list.reaction_card(reaction_name)).to_be_visible()


def test_reaction_names_must_be_unique(logged_in_page: Page) -> None:
    """A name already used in the workbook is refused, and the user stays on the list."""
    reaction_list = ReactionListPage(logged_in_page)
    reaction_list.goto(TEST_WORKGROUP, TEST_WORKBOOK)

    reaction_list.create_reaction(SEEDED_REACTION_NAME)

    expect(reaction_list.new_reaction_error).to_contain_text("already")
    expect(logged_in_page).not_to_have_url(re.compile(r"/sketcher/"))


def test_sketched_reaction_fills_and_keeps_reaction_table(
    logged_in_page: Page,
) -> None:
    """
    Drawing a reaction lists its compounds in the reaction table. Once the
    table has been edited and autosaved, reopening the reaction from its
    workbook brings back the same table.
    """
    reaction_name = "E2E Sketched Reaction"
    description = "Phosphination of the test family"

    reaction_list = ReactionListPage(logged_in_page)
    reaction_list.goto(TEST_WORKGROUP, TEST_WORKBOOK)
    reaction_list.create_reaction(reaction_name)

    reaction_page = ReactionConstructorPage(logged_in_page)
    expect(reaction_page.heading).to_contain_text(reaction_name)
    reaction_page.draw_reaction(REACTION_SMILES)

    reaction_page.wait_for_reaction_table()
    for number, reactant in enumerate(REACTANT_NAMES, start=1):
        expect(reaction_page.reactant_name(number)).to_have_value(reactant)
    expect(reaction_page.product_name(1)).to_have_value(PRODUCT_NAME)
    reaction_page.describe_reaction(description)

    reaction_list.goto(TEST_WORKGROUP, TEST_WORKBOOK)
    reaction_list.reload_reaction(reaction_name)

    expect(reaction_page.heading).to_contain_text(reaction_name)
    reaction_page.wait_for_reaction_table()
    expect(reaction_page.reaction_description).to_have_value(description)
    for number, reactant in enumerate(REACTANT_NAMES, start=1):
        expect(reaction_page.reactant_name(number)).to_have_value(reactant)
    expect(reaction_page.product_name(1)).to_have_value(PRODUCT_NAME)


def test_user_deletes_own_reaction(logged_in_page: Page) -> None:
    """After confirming, a deleted reaction is no longer listed in its workbook."""
    reaction_name = "E2E Reaction To Delete"

    reaction_list = ReactionListPage(logged_in_page)
    reaction_list.goto(TEST_WORKGROUP, TEST_WORKBOOK)
    reaction_list.create_reaction(reaction_name)
    expect(ReactionConstructorPage(logged_in_page).heading).to_contain_text(
        reaction_name
    )

    reaction_list.goto(TEST_WORKGROUP, TEST_WORKBOOK)
    reaction_list.delete_reaction(reaction_name)

    expect(reaction_list.reaction_card(reaction_name)).to_have_count(0)
    expect(reaction_list.reaction_card(SEEDED_REACTION_NAME)).to_be_visible()
