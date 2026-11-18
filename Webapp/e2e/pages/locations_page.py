"""Page object for a workgroup's inventory storage locations in the React app."""

from e2e.pages.base_page import BasePage
from playwright.sync_api import Locator


class LocationsPage(BasePage):
    """The locations list at /spa/inventory/<workgroup>/locations."""

    @property
    def heading(self) -> Locator:
        """The "Locations" page heading."""
        return self.page.get_by_role("heading", name="Locations", exact=True)

    @property
    def add_location_dialog(self) -> Locator:
        """The "Add location" modal, while it is open."""
        return self.page.get_by_role("dialog", name="Add location")

    def add_location(self, name: str, description: str = "") -> None:
        """
        Open "+ Add location", fill in the form and submit it.

        Args:
            name: Location name, unique within the workgroup.
            description: Optional description of where it is.
        """
        self.page.get_by_role("button", name="+ Add location").click()
        dialog = self.add_location_dialog
        dialog.get_by_label("Name").fill(name)
        if description:
            dialog.get_by_label("Description").fill(description)
        dialog.get_by_role("button", name="Add", exact=True).click()

    def location_row(self, name: str) -> Locator:
        """
        Find a location's row.

        Args:
            name: The location's name.

        Returns:
            Locator for the row, which holds its container count and status.
        """
        return self.page.get_by_role("row").filter(
            has=self.page.get_by_role("button", name=f"Edit {name}", exact=True)
        )
