"""
Journeys: H&S admins design COSHH templates, and users fill in COSHH forms
made from them, standalone or for a reaction.

The template, its section and its text field are seeded by
e2e/populate_database.py, because building one in the form designer is a
separate journey from filling one in.
"""

import re
from typing import Callable

from e2e.pages.coshh_pages import CoshhFormManagementPage, MyCoshhFormsPage
from e2e.pages.reaction_constructor_page import ReactionConstructorPage
from e2e.pages.reaction_list_page import ReactionListPage
from e2e.seed_data import (
    COSHH_SECTION,
    COSHH_TEMPLATE,
    COSHH_TEXT_FIELD,
    HS_ADMIN_USERNAME,
    ROLE_USER_PASSWORD,
    TEST_WORKBOOK,
    TEST_WORKGROUP,
)
from playwright.sync_api import Page, expect


def test_hs_admin_creates_template(log_in_as: Callable[[str, str], Page]) -> None:
    """An H&S admin creates a COSHH template, which is then listed as a draft."""
    template_name = "E2E Admin Template"
    page = log_in_as(HS_ADMIN_USERNAME, ROLE_USER_PASSWORD)

    management = CoshhFormManagementPage(page)
    management.goto()
    management.create_template(template_name, description="For solvent handling")

    expect(management.toast("Project Type created")).to_be_visible()
    row = management.template_row(template_name)
    expect(row).to_contain_text("For solvent handling")
    expect(row).to_contain_text("Draft")


def test_template_needs_a_description(log_in_as: Callable[[str, str], Page]) -> None:
    """The create template form won't submit without a description."""
    page = log_in_as(HS_ADMIN_USERNAME, ROLE_USER_PASSWORD)

    management = CoshhFormManagementPage(page)
    management.goto()
    management.create_template("E2E Undescribed Template", description="")

    expect(management.create_dialog).to_contain_text(
        "Project type description required"
    )


def test_standard_user_cannot_create_templates(logged_in_page: Page) -> None:
    """Users who aren't H&S admins can see the templates but not create one."""
    management = CoshhFormManagementPage(logged_in_page)
    management.goto()

    expect(management.template_row(COSHH_TEMPLATE)).to_be_visible()
    expect(management.new_button).to_have_count(0)


def test_user_fills_in_standalone_coshh_form(logged_in_page: Page) -> None:
    """
    A user starts a COSHH form from My COSHH Forms. What they type is saved
    automatically, kept after a reload, and the form is listed as theirs.
    """
    control_measures = "Use in a fume hood with nitrile gloves"

    my_forms = MyCoshhFormsPage(logged_in_page)
    my_forms.goto()
    my_forms.create_form(COSHH_TEMPLATE)

    expect(my_forms.toast("COSHH form created")).to_be_visible()
    expect(logged_in_page).to_have_url(re.compile(r"/spa/my-coshh-forms/\d+$"))
    form = my_forms.form
    expect(form.section(COSHH_SECTION)).to_be_visible()

    form.fill_and_wait_for_autosave(COSHH_TEXT_FIELD, control_measures)
    logged_in_page.reload()
    expect(form.text_field(COSHH_TEXT_FIELD)).to_have_value(control_measures)

    my_forms.goto()
    expect(my_forms.form_rows(COSHH_TEMPLATE).first).to_be_visible()


def test_user_adds_coshh_form_to_reaction(logged_in_page: Page) -> None:
    """
    A user creates a COSHH form from a reaction's page. It is saved against the
    reaction, so the form is still there when the reaction is opened again.
    """
    reaction_name = "E2E Reaction With COSHH"
    control_measures = "Work behind a blast shield"

    reaction_list = ReactionListPage(logged_in_page)
    reaction_list.goto(TEST_WORKGROUP, TEST_WORKBOOK)
    reaction_list.create_reaction(reaction_name)

    reaction_page = ReactionConstructorPage(logged_in_page)
    expect(reaction_page.heading).to_contain_text(reaction_name)
    reaction_page.open_coshh_section()
    reaction_page.create_coshh_form(COSHH_TEMPLATE)

    expect(reaction_page.toast("COSHH form created")).to_be_visible()
    form = reaction_page.coshh_form
    expect(form.section(COSHH_SECTION)).to_be_visible()
    form.fill_and_wait_for_autosave(COSHH_TEXT_FIELD, control_measures)

    logged_in_page.reload()
    reaction_page.open_coshh_section()
    expect(form.text_field(COSHH_TEXT_FIELD)).to_have_value(control_measures)
