"""Page object for the reactions listed in a workbook on the workgroup page."""

from urllib.parse import quote

from e2e.pages.base_page import BasePage
from playwright.sync_api import Locator


class ReactionListPage(BasePage):
    """The workgroup page at /workgroup/<workgroup>/<workbook>, with the workbook's reaction cards."""

    def goto(self, workgroup: str, workbook: str) -> None:
        """
        Open the workgroup page with a workbook selected.

        Args:
            workgroup: Name of the workgroup.
            workbook: Name of the workbook whose reactions to list.
        """
        self.page.goto(f"/workgroup/{quote(workgroup)}/{quote(workbook)}")
        # The reaction list is fetched after the page loads, and brings a second
        # copy of the "New Reaction" modal. Until it arrives, a name typed into
        # the modal is read from the wrong copy, so wait for the list's hidden
        # reaction count, which comes with it.
        self.page.locator("#js-number-of-reactions").wait_for(state="attached")

    @property
    def new_reaction_dialog(self) -> Locator:
        """The "New Reaction" modal, while it is open."""
        return self.page.get_by_role("dialog", name="New Reaction")

    @property
    def new_reaction_error(self) -> Locator:
        """
        The error message in the "New Reaction" modal.

        It has no role or label in the template, so it is found by id. The
        page renders the modal twice, so only the open one is searched.
        """
        return self.new_reaction_dialog.locator("#error-warning-new-reaction")

    def create_reaction(self, name: str) -> None:
        """
        Open "+ New Reaction", name the reaction and submit.

        On success the browser goes to the new reaction's page.

        Args:
            name: Reaction name, unique within the workbook.
        """
        self.page.get_by_role("button", name="+ New Reaction").click()
        dialog = self.new_reaction_dialog
        dialog.get_by_label("Reaction Name").fill(name)
        dialog.get_by_role("button", name="Create").click()

    def reaction_card(self, name: str) -> Locator:
        """
        Find a reaction's card.

        Args:
            name: The reaction's name.

        Returns:
            Locator for the card, which holds its ID and Reload/Delete/Clone buttons.
        """
        return self.page.locator(".reaction-card").filter(
            has=self.page.get_by_role("heading", name=name, exact=True)
        )

    def reload_reaction(self, name: str) -> None:
        """
        Click a reaction's "Reload" button to open it.

        Args:
            name: The reaction's name.
        """
        self.reaction_card(name).get_by_role("button", name="Reload").click()

    def delete_reaction(self, name: str) -> None:
        """
        Click a reaction's "Delete" button and accept the browser's confirmation.

        Args:
            name: The reaction's name.
        """
        self.page.once("dialog", lambda dialog: dialog.accept())
        self.reaction_card(name).get_by_role("button", name="Delete").click()
