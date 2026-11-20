"""
Journeys: a user records a reaction's workup and purification step by step on
the reaction page, finds them as they left them when the reaction is reopened,
saves the steps as a workbook template, and starts another reaction's workup
from that template.
"""

from e2e.pages.reaction_constructor_page import ReactionConstructorPage
from e2e.pages.reaction_list_page import ReactionListPage
from e2e.pages.workup_section import WorkupSection
from e2e.seed_data import TEST_WORKBOOK, TEST_WORKGROUP
from playwright.sync_api import Page, expect


def _create_reaction(page: Page, name: str) -> WorkupSection:
    """Creates a reaction, opens its notes and returns its workup card."""
    reaction_list = ReactionListPage(page)
    reaction_list.goto(TEST_WORKGROUP, TEST_WORKBOOK)
    reaction_list.create_reaction(name)
    expect(ReactionConstructorPage(page).heading).to_contain_text(name)

    workup = WorkupSection(page)
    workup.open()
    return workup


def test_user_records_a_wash_step(logged_in_page: Page) -> None:
    """
    A user adds a wash to a reaction's workup, which opens ready to fill in, picks
    the organic phase from the coloured CHEM21 solvent list and types an aqueous
    phase that isn't in the database. The step saves itself, so both are still
    there, and the solvent still linked, when the reaction is opened again.
    """
    workup = _create_reaction(logged_in_page, "E2E Reaction With Workup")

    wash = workup.add_step("Wash")
    expect(wash.get_by_role("button", name="Collapse step")).to_be_visible()
    # Testanol is a seeded CHEM21 solvent flagged "Hazardous"
    workup.choose_solvent_and_wait_for_save(wash, "Organic phase", "Testa", "Testanol")
    workup.fill_and_wait_for_save(wash, "Aqueous phase", "Brine")

    logged_in_page.reload()
    workup.open()
    wash = workup.step(1, "Wash")
    workup.expand_step(wash)
    expect(wash.get_by_label("Organic phase", exact=True)).to_have_value("Testanol")
    expect(wash.get_by_text("Hazardous", exact=True)).to_be_visible()
    expect(wash.get_by_label("Aqueous phase", exact=True)).to_have_value("Brine")
    expect(wash.get_by_text("Not in database", exact=True)).to_be_visible()


def test_saved_workup_renders_when_reaction_is_reopened(logged_in_page: Page) -> None:
    """
    A user records two workup steps and a purification step, then reopens the
    reaction from its workbook. Each card shows its own steps in the order they
    were added, collapsed to a summary of what was entered, and expanding a step
    shows every saved value back in its input.
    """
    reaction_name = "E2E Reopened Workup"
    workup = _create_reaction(logged_in_page, reaction_name)
    purification = WorkupSection(logged_in_page, "PURIFICATION")

    wash = workup.add_step("Wash")
    workup.choose_solvent_and_wait_for_save(wash, "Organic phase", "Testa", "Testanol")
    workup.fill_and_wait_for_save(wash, "Aqueous phase", "Brine")
    workup.fill_and_wait_for_save(wash, "Number of washes", "3")

    concentration = workup.add_step("Concentration")
    workup.select_and_wait_for_save(concentration, "Method", "Rotary evaporator")
    workup.fill_and_wait_for_save(concentration, "Bath temperature", "40")
    workup.fill_and_wait_for_save(concentration, "Pressure", "50")

    column = purification.add_step("Flash chromatography")
    purification.fill_and_wait_for_save(column, "Gradient", "0-20% EtOAc")

    reaction_list = ReactionListPage(logged_in_page)
    reaction_list.goto(TEST_WORKGROUP, TEST_WORKBOOK)
    reaction_list.reload_reaction(reaction_name)
    expect(ReactionConstructorPage(logged_in_page).heading).to_contain_text(
        reaction_name
    )
    workup.open()

    # collapsed, each step shows a summary of what was entered
    expect(workup.steps).to_have_count(2)
    wash = workup.step(1, "Wash")
    concentration = workup.step(2, "Concentration")
    for step in (wash, concentration):
        expect(step.get_by_role("button", name="Expand step")).to_be_visible()
    expect(wash).to_contain_text("Testanol / 3 × Brine")
    expect(concentration).to_contain_text("Rotary evaporator · 40 °C · 50 mbar")
    expect(purification.steps).to_have_count(1)
    column = purification.step(1, "Flash chromatography")
    expect(column).to_contain_text("0-20% EtOAc")

    workup.expand_step(wash)
    expect(wash.get_by_label("Organic phase", exact=True)).to_have_value("Testanol")
    expect(wash.get_by_text("Hazardous", exact=True)).to_be_visible()
    expect(wash.get_by_label("Aqueous phase", exact=True)).to_have_value("Brine")
    expect(wash.get_by_label("Number of washes", exact=True)).to_have_value("3")

    workup.expand_step(concentration)
    expect(concentration.get_by_label("Method", exact=True)).to_have_value(
        "Rotary evaporator"
    )
    expect(concentration.get_by_label("Bath temperature", exact=True)).to_have_value(
        "40"
    )
    expect(concentration.get_by_label("Pressure", exact=True)).to_have_value("50")

    purification.expand_step(column)
    expect(column.get_by_label("Gradient", exact=True)).to_have_value("0-20% EtOAc")


def test_user_reuses_a_workup_as_a_template(logged_in_page: Page) -> None:
    """
    A user saves a reaction's workup steps as a template and starts another
    reaction's workup from it. The new reaction gets the steps without the values.
    """
    template_name = "E2E Aqueous Workup"
    workup = _create_reaction(logged_in_page, "E2E Workup Template Source")
    wash = workup.add_step("Wash")
    workup.fill_and_wait_for_save(wash, "Organic phase", "DCM")
    expect(workup.add_step("Dry with material")).to_be_visible()

    workup.save_as_template(template_name)
    expect(workup.toast("Template saved")).to_be_visible()

    workup = _create_reaction(logged_in_page, "E2E Workup Template Target")
    workup.start_from_template(template_name)

    expect(workup.toast(f'Added the steps from "{template_name}"')).to_be_visible()
    wash = workup.step(1, "Wash")
    expect(wash.get_by_label("Organic phase", exact=True)).to_have_value("")
    expect(workup.step(2, "Dry with material")).to_be_visible()
