"""Page object for the create workbook form."""

from e2e.pages.base_page import BasePage


class CreateWorkbookPage(BasePage):
    """The form at /create_workbook/<workgroup>."""

    def goto(self, workgroup: str) -> None:
        """
        Open the form directly by URL, as a user could without using the UI.

        Args:
            workgroup: Name of the workgroup to create the workbook in.
        """
        self.page.goto(f"/create_workbook/{workgroup}")

    def create(self, name: str, abbreviation: str) -> None:
        """
        Fill in and submit the form.

        The inputs have no <label> elements in the template, so they are found by id.

        Args:
            name: Workbook name, at least 4 characters.
            abbreviation: 3 character alphanumeric abbreviation.
        """
        self.page.locator("#workbook").fill(name)
        self.page.locator("#abbreviation").fill(abbreviation)
        self.page.get_by_role("button", name="Create Workbook").click()
