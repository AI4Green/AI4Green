"""Page object for the logged-in home page and its Quick Access panel."""

import re

from e2e.pages.base_page import BasePage
from playwright.sync_api import Locator


class HomePage(BasePage):
    """The home page, where users pick a workgroup and workbook from Quick Access."""

    PATH = "/home"

    def goto(self) -> None:
        """Open the home page."""
        self.page.goto(self.PATH)

    @property
    def welcome_heading(self) -> Locator:
        """The "Welcome to AI4Green, <username>!" heading."""
        return self.page.get_by_role("heading", name=re.compile("Welcome to AI4Green"))

    def select_workgroup(self, workgroup_name: str) -> None:
        """
        Click a workgroup in Quick Access, which loads its workbooks.

        Args:
            workgroup_name: Name of the workgroup.
        """
        self.page.get_by_role("button", name=workgroup_name, exact=True).click()

    def workbook_button(self, workbook_name: str) -> Locator:
        """
        Find a workbook in Quick Access. Only shown once a workgroup is selected.

        Args:
            workbook_name: Name of the workbook.

        Returns:
            Locator for the workbook's button.
        """
        return self.page.get_by_role("button", name=workbook_name, exact=True)

    def start_new_workbook(self) -> None:
        """Click "New Workbook" in Quick Access. Only shown once a workgroup is selected."""
        self.page.get_by_role("button", name="New Workbook", exact=True).click()
