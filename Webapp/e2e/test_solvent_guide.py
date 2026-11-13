"""Journey: a user looks up a solvent in the CHEM21 solvent guide."""

from e2e.pages.solvent_guide_page import SolventGuidePage
from playwright.sync_api import Page, expect


def test_user_views_solvent_card(logged_in_page: Page) -> None:
    """Choosing a solvent shows its card with its CHEM21 ranking."""
    solvent_guide = SolventGuidePage(logged_in_page)
    solvent_guide.goto()

    solvent_guide.choose_solvent(family="Alcohols", solvent="Ethanol")

    card = solvent_guide.solvent_card("Ethanol")
    expect(card).to_be_visible()
    expect(card).to_contain_text("Recommended")
