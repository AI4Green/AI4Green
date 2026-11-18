"""
Journeys: managing a workgroup's chemical inventory, and who is allowed to.

The inventory workgroup is seeded with two locations and one stocked compound
(see e2e/populate_database.py). Tests that act on a container get a fresh one
from the ``inventory_container`` fixture, so they don't depend on each other.
"""

from typing import Callable

from e2e.pages.inventory_page import InventoryPage
from e2e.pages.locations_page import LocationsPage
from e2e.seed_data import (
    INVENTORY_LOCATION,
    INVENTORY_OTHER_LOCATION,
    INVENTORY_WORKGROUP,
    ROLE_USER_PASSWORD,
    STANDARD_MEMBER_USERNAME,
    STOCKED_COMPOUND_NAME,
    UNSTOCKED_COMPOUND_CAS,
    UNSTOCKED_COMPOUND_NAME,
)
from playwright.sync_api import Page, expect


def test_inventory_opens_from_workgroup_page(logged_in_page: Page) -> None:
    """The Inventory button on the workgroup page opens the React inventory."""
    logged_in_page.goto(f"/workgroup/{INVENTORY_WORKGROUP}")
    logged_in_page.get_by_role("link", name="Inventory", exact=True).click()

    inventory = InventoryPage(logged_in_page)
    expect(inventory.heading).to_be_visible()
    expect(inventory.compound_row(STOCKED_COMPOUND_NAME)).to_be_visible()


def test_pi_adds_location(logged_in_page: Page) -> None:
    """The PI adds a storage location, which is then listed as active and empty."""
    location = "E2E Fridge"

    inventory = InventoryPage(logged_in_page)
    inventory.goto(INVENTORY_WORKGROUP)
    inventory.manage_locations()

    locations = LocationsPage(logged_in_page)
    expect(locations.heading).to_be_visible()
    locations.add_location(location, description="Room B12, bottom shelf")

    expect(locations.toast("Location added")).to_be_visible()
    row = locations.location_row(location)
    expect(row).to_contain_text("Room B12, bottom shelf")
    expect(row).to_contain_text("Active")


def test_location_names_must_be_unique(logged_in_page: Page) -> None:
    """Reusing a location's name, in any case, is refused before it reaches the server."""
    InventoryPage(logged_in_page).goto(INVENTORY_WORKGROUP)
    InventoryPage(logged_in_page).manage_locations()

    locations = LocationsPage(logged_in_page)
    locations.add_location(INVENTORY_LOCATION.upper())

    expect(locations.add_location_dialog).to_contain_text(
        "There's already a location with this name"
    )


def test_pi_adds_stock_of_new_compound(logged_in_page: Page) -> None:
    """
    The PI finds a compound the workgroup doesn't stock yet by its CAS number
    and adds a container of it.
    """
    lot_number = "E2E-NEW-STOCK"

    inventory = InventoryPage(logged_in_page)
    inventory.goto(INVENTORY_WORKGROUP)
    expect(inventory.compound_row(UNSTOCKED_COMPOUND_NAME)).to_have_count(0)

    inventory.add_stock(
        cas=UNSTOCKED_COMPOUND_CAS,
        location=INVENTORY_LOCATION,
        quantity="25",
        unit="g",
        lot_number=lot_number,
    )

    expect(inventory.toast("Stock added")).to_be_visible()
    row = inventory.compound_row(UNSTOCKED_COMPOUND_NAME)
    expect(row).to_contain_text("25 g")
    expect(row).to_contain_text(INVENTORY_LOCATION)

    inventory.expand_compound(UNSTOCKED_COMPOUND_NAME)
    expect(inventory.container_row(lot_number)).to_contain_text(INVENTORY_LOCATION)


def test_pi_adjusts_container_quantity(
    logged_in_page: Page, inventory_container: str
) -> None:
    """After a stock check, the PI sets a container's quantity to what was measured."""
    inventory = InventoryPage(logged_in_page)
    inventory.goto(INVENTORY_WORKGROUP)
    inventory.expand_compound(STOCKED_COMPOUND_NAME)
    expect(inventory.container_row(inventory_container)).to_contain_text("500 mL")

    inventory.adjust_container(inventory_container, quantity="320", notes="Stock check")

    expect(inventory.toast("Quantity updated")).to_be_visible()
    expect(inventory.container_row(inventory_container)).to_contain_text("320 mL")


def test_pi_moves_container(logged_in_page: Page, inventory_container: str) -> None:
    """The PI moves a container to the workgroup's other location."""
    inventory = InventoryPage(logged_in_page)
    inventory.goto(INVENTORY_WORKGROUP)
    inventory.expand_compound(STOCKED_COMPOUND_NAME)

    inventory.move_container(inventory_container, INVENTORY_OTHER_LOCATION)

    expect(inventory.toast("Container moved")).to_be_visible()
    expect(inventory.container_row(inventory_container)).to_contain_text(
        INVENTORY_OTHER_LOCATION
    )


def test_pi_disposes_of_container(
    logged_in_page: Page, inventory_container: str
) -> None:
    """
    A disposed container drops out of the list, and comes back as archived
    when "Show archived" is ticked.
    """
    inventory = InventoryPage(logged_in_page)
    inventory.goto(INVENTORY_WORKGROUP)
    inventory.expand_compound(STOCKED_COMPOUND_NAME)

    inventory.dispose_of_container(inventory_container, reason="Contaminated")

    expect(inventory.toast("Container disposed")).to_be_visible()
    expect(inventory.container_row(inventory_container)).to_have_count(0)

    inventory.show_archived()
    archived_row = inventory.container_row(inventory_container)
    expect(archived_row).to_contain_text("Archived")
    # an archived container can only show its history
    expect(inventory.container_action("Adjust", inventory_container)).to_have_count(0)
    expect(inventory.container_action("History", inventory_container)).to_be_visible()


def test_standard_member_can_view_but_not_change_stock(
    log_in_as: Callable[[str, str], Page], inventory_container: str
) -> None:
    """A standard member sees the stock and its history, but none of the actions that change it."""
    page = log_in_as(STANDARD_MEMBER_USERNAME, ROLE_USER_PASSWORD)

    inventory = InventoryPage(page)
    inventory.goto(INVENTORY_WORKGROUP)
    inventory.expand_compound(STOCKED_COMPOUND_NAME)

    expect(inventory.container_action("History", inventory_container)).to_be_visible()
    expect(inventory.container_action("Adjust", inventory_container)).to_have_count(0)
    expect(inventory.container_action("Dispose", inventory_container)).to_have_count(0)
    expect(inventory.add_stock_button).to_have_count(0)
