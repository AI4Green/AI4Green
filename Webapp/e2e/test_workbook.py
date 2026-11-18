"""Journey: creating a workbook from the home page, and who is allowed to."""

from typing import Callable

from e2e.pages.create_workbook_page import CreateWorkbookPage
from e2e.pages.home_page import HomePage
from e2e.seed_data import (
    ROLE_USER_PASSWORD,
    ROLES_WORKGROUP,
    SENIOR_RESEARCHER_USERNAME,
    STANDARD_MEMBER_USERNAME,
    TEST_WORKBOOK,
    TEST_WORKGROUP,
)
from playwright.sync_api import Page, expect

PERMISSION_DENIED = "You do not have permission to view this page"


def test_seeded_workbook_shows_in_quick_access(logged_in_page: Page) -> None:
    """Selecting a workgroup on the home page lists its workbooks."""
    home_page = HomePage(logged_in_page)
    home_page.goto()

    home_page.select_workgroup(TEST_WORKGROUP)

    expect(home_page.workbook_button(TEST_WORKBOOK)).to_be_visible()


def test_pi_creates_workbook(logged_in_page: Page) -> None:
    """The PI creates a workbook from Quick Access and then sees it listed."""
    new_workbook = "E2E Workbook"

    home_page = HomePage(logged_in_page)
    home_page.goto()
    home_page.select_workgroup(TEST_WORKGROUP)
    home_page.start_new_workbook()

    create_page = CreateWorkbookPage(logged_in_page)
    create_page.create(name=new_workbook, abbreviation="E2E")

    expect(create_page.flash_message("Workbook has been created!")).to_be_visible()

    home_page.goto()
    home_page.select_workgroup(TEST_WORKGROUP)
    expect(home_page.workbook_button(new_workbook)).to_be_visible()


def test_senior_researcher_creates_workbook(
    log_in_as: Callable[[str, str], Page]
) -> None:
    """A senior researcher is allowed to create a workbook in their workgroup."""
    new_workbook = "Senior Workbook"
    page = log_in_as(SENIOR_RESEARCHER_USERNAME, ROLE_USER_PASSWORD)

    home_page = HomePage(page)
    home_page.select_workgroup(ROLES_WORKGROUP)
    home_page.start_new_workbook()

    create_page = CreateWorkbookPage(page)
    create_page.create(name=new_workbook, abbreviation="SRW")

    expect(create_page.flash_message("Workbook has been created!")).to_be_visible()


def test_standard_member_cannot_create_workbook(
    log_in_as: Callable[[str, str], Page]
) -> None:
    """
    A standard member who goes straight to the create workbook URL is sent home
    with a permission message, and never sees the form.
    """
    page = log_in_as(STANDARD_MEMBER_USERNAME, ROLE_USER_PASSWORD)

    create_page = CreateWorkbookPage(page)
    create_page.goto(ROLES_WORKGROUP)

    expect(create_page.flash_message(PERMISSION_DENIED)).to_be_visible()
    expect(HomePage(page).welcome_heading).to_be_visible()
    expect(page.get_by_role("button", name="Create Workbook")).to_be_hidden()
