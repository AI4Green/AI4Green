"""Page object for the sign-in page."""

from e2e.pages.base_page import BasePage


class LoginPage(BasePage):
    """The username/password sign-in form at /auth/login."""

    PATH = "/auth/login"

    def goto(self) -> None:
        """Open the sign-in page."""
        self.page.goto(self.PATH)

    def login(self, username: str, password: str) -> None:
        """
        Fill in and submit the sign-in form.

        Args:
            username: Username or email address.
            password: Password.
        """
        self.page.get_by_label("Username or Email").fill(username)
        self.page.get_by_label("Password", exact=True).fill(password)
        self.page.get_by_role("button", name="Log In").click()
