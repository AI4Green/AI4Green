"""Page object for a workgroup's inventory in the React app."""

from urllib.parse import quote

from e2e.pages.base_page import BasePage
from playwright.sync_api import Locator


class InventoryPage(BasePage):
    """The inventory at /spa/inventory/<workgroup>: compounds, their containers and stock actions."""

    def goto(self, workgroup: str) -> None:
        """
        Open a workgroup's inventory.

        Args:
            workgroup: Name of the workgroup.
        """
        self.page.goto(f"/spa/inventory/{quote(workgroup)}")

    @property
    def heading(self) -> Locator:
        """The "Inventory" page heading."""
        return self.page.get_by_role("heading", name="Inventory", exact=True)

    @property
    def add_stock_button(self) -> Locator:
        """The "+ Add stock" button, shown only to users who can change stock."""
        return self.page.get_by_role("button", name="+ Add stock")

    def dialog(self, title: str) -> Locator:
        """
        Find an open modal by its title.

        Args:
            title: The modal's heading, e.g. "Adjust quantity".

        Returns:
            Locator for the modal.
        """
        return self.page.get_by_role("dialog", name=title)

    def add_stock(
        self, cas: str, location: str, quantity: str, unit: str, lot_number: str
    ) -> None:
        """
        Add a container of a compound found by CAS number with "+ Add stock".

        Args:
            cas: CAS number to search the compound database for.
            location: Name of the location to store the container in.
            quantity: Amount in the container.
            unit: Unit of the quantity, e.g. "g".
            lot_number: Lot number printed on the container.
        """
        self.add_stock_button.click()
        dialog = self.dialog("Add stock")
        dialog.get_by_placeholder("Name or CAS number").fill(cas)
        dialog.get_by_role("button", name="Search database").click()
        dialog.get_by_role("button", name="Use this compound").click()
        dialog.get_by_label("Location").select_option(label=location)
        dialog.get_by_label("Quantity").fill(quantity)
        dialog.get_by_label("Unit").select_option(label=unit)
        dialog.get_by_label("Lot number").fill(lot_number)
        dialog.get_by_role("button", name="Add", exact=True).click()

    def compound_row(self, compound_name: str) -> Locator:
        """
        Find a compound's row in the compound view.

        Args:
            compound_name: The compound's name.

        Returns:
            Locator for the row, which holds its stock totals and locations.
        """
        return self.page.get_by_role("row").filter(
            has=self.page.get_by_role("button", name=f"Expand {compound_name}")
        )

    def expand_compound(self, compound_name: str) -> None:
        """
        Expand a compound's row to list its containers.

        Args:
            compound_name: The compound's name.
        """
        self.page.get_by_role("button", name=f"Expand {compound_name}").click()

    def container_row(self, lot_number: str) -> Locator:
        """
        Find a container's row. Its compound must be expanded first.

        Args:
            lot_number: The container's lot number.

        Returns:
            Locator for the row.
        """
        # an expanded compound's containers are a table nested in one of the
        # compound table's rows, so skip rows that hold other rows
        return self.page.get_by_role("row").filter(
            has=self.page.get_by_role("cell", name=lot_number, exact=True),
            has_not=self.page.get_by_role("row"),
        )

    def container_action(self, action: str, lot_number: str) -> Locator:
        """
        Find one of a container's action buttons. Its compound must be expanded first.

        Args:
            action: "Adjust", "Move", "Dispose" or "History".
            lot_number: The container's lot number.

        Returns:
            Locator for the button.
        """
        return self.page.get_by_role(
            "button", name=f"{action} container {lot_number}", exact=True
        )

    def adjust_container(self, lot_number: str, quantity: str, notes: str) -> None:
        """
        Set a container's quantity to a measured amount.

        Args:
            lot_number: The container's lot number.
            quantity: The new quantity, in the container's unit.
            notes: Why the quantity changed.
        """
        self.container_action("Adjust", lot_number).click()
        dialog = self.dialog("Adjust quantity")
        dialog.get_by_label("New quantity").fill(quantity)
        dialog.get_by_label("Notes").fill(notes)
        dialog.get_by_role("button", name="Save").click()

    def move_container(self, lot_number: str, location: str) -> None:
        """
        Move a container to another location.

        Args:
            lot_number: The container's lot number.
            location: Name of the location to move it to.
        """
        self.container_action("Move", lot_number).click()
        dialog = self.dialog("Move container")
        dialog.get_by_label("New location").select_option(label=location)
        dialog.get_by_role("button", name="Move", exact=True).click()

    def dispose_of_container(self, lot_number: str, reason: str) -> None:
        """
        Dispose of a container, which empties and archives it.

        Args:
            lot_number: The container's lot number.
            reason: Why it was disposed of.
        """
        self.container_action("Dispose", lot_number).click()
        dialog = self.dialog("Dispose of container")
        dialog.get_by_label("Reason").fill(reason)
        dialog.get_by_role("button", name="Dispose", exact=True).click()

    def show_archived(self) -> None:
        """Tick "Show archived" so archived containers are listed too."""
        self.page.get_by_text("Show archived", exact=True).click()

    def manage_locations(self) -> None:
        """Click "Manage locations" to go to the locations page."""
        self.page.get_by_role("link", name="Manage locations").click()
