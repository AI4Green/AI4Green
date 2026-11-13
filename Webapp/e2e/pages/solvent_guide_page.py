"""Page object for the CHEM21 solvent guide."""

from e2e.pages.base_page import BasePage
from playwright.sync_api import Locator


class SolventGuidePage(BasePage):
    """The solvent guide, where solvents are grouped by family and shown as cards."""

    PATH = "/solvent_guide"

    def goto(self) -> None:
        """Open the solvent guide."""
        self.page.goto(self.PATH)

    def choose_solvent(self, family: str, solvent: str) -> None:
        """
        Expand a solvent family and click one of its solvents.

        Args:
            family: Family name, e.g. "Alcohols".
            solvent: Solvent name, e.g. "Ethanol".
        """
        self.page.get_by_role("button", name=family, exact=True).click()
        self.page.get_by_role("button", name=solvent, exact=True).click()

    def solvent_card(self, solvent: str) -> Locator:
        """
        Find the visible card showing a solvent's details.

        Args:
            solvent: Solvent name shown as the card title.

        Returns:
            Locator for the visible tab panel holding that card.
        """
        return self.page.get_by_role("tabpanel").filter(
            has=self.page.get_by_role("heading", name=solvent, exact=True)
        )
