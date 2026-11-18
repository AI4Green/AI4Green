"""Journey: a user signs in and out."""

from e2e.pages.home_page import HomePage
from e2e.pages.login_page import LoginPage
from e2e.seed_data import TEST_PASSWORD, TEST_USERNAME
from playwright.sync_api import Page, expect


def test_user_can_log_in_and_out(page: Page) -> None:
    """A verified user signs in, sees the home page, then logs out."""
    login_page = LoginPage(page)
    login_page.goto()
    login_page.login(TEST_USERNAME, TEST_PASSWORD)

    home_page = HomePage(page)
    expect(home_page.welcome_heading).to_contain_text(TEST_USERNAME)

    home_page.logout()

    expect(home_page.welcome_heading).to_be_hidden()
    expect(home_page.login_link).to_be_visible()


def test_wrong_password_shows_error(page: Page) -> None:
    """A wrong password keeps the user logged out and shows a flash message."""
    login_page = LoginPage(page)
    login_page.goto()
    login_page.login(TEST_USERNAME, "not-the-password")

    expect(login_page.flash_message("Invalid username or password")).to_be_visible()
    expect(login_page.login_link).to_be_visible()
