"""Parts of the page that every AI4Green page shares: the navbar and flash messages."""

from playwright.sync_api import Locator, Page


class BasePage:
    """Shared navbar and flash message behaviour."""

    def __init__(self, page: Page) -> None:
        """
        Args:
            page: The Playwright page to drive.
        """
        self.page = page

    def flash_message(self, text: str) -> Locator:
        """
        Find a flashed alert containing the given text.

        Args:
            text: Text the alert should contain.

        Returns:
            Locator for the matching alert.
        """
        return self.page.get_by_role("alert").filter(has_text=text)

    @property
    def login_link(self) -> Locator:
        """The navbar "Log in" link, shown only to logged-out users."""
        return self.page.get_by_role("link", name="Log in")

    def logout(self) -> None:
        """Open the user menu in the navbar and click Logout."""
        self.page.locator("#user-dropdown").click()
        self.page.get_by_role("link", name="Logout").click()
