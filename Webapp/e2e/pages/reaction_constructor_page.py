"""Page object for a reaction's page: the sketcher, reaction table and COSHH form."""

from e2e.pages.base_page import BasePage
from e2e.pages.coshh_pages import CoshhForm, choose_template
from playwright.sync_api import Locator
from playwright.sync_api import TimeoutError as PlaywrightTimeoutError

# Ketcher is a large app loaded into an iframe, so allow longer than usual
KETCHER_LOAD_TIMEOUT_MS = 30_000
# the page autosaves a sketch 800 ms after the last change
TABLE_REQUEST_TIMEOUT_MS = 5_000
DRAW_ATTEMPTS = 3
REACTION_TABLE_TIMEOUT_MS = 30_000


class ReactionConstructorPage(BasePage):
    """
    The reaction page at /sketcher/<workgroup>/<workbook>/<reaction id>/no.

    Drawing the reaction updates the reaction table: the page autosaves the
    sketch and the server renders the table for the compounds it recognises.
    """

    @property
    def heading(self) -> Locator:
        """The page heading, which shows the reaction's ID and name."""
        return self.page.get_by_role("heading", level=2)

    def draw_reaction(self, reaction_smiles: str) -> None:
        """
        Put a reaction in the Ketcher sketcher.

        Drawing on Ketcher's canvas with the mouse is impractical to automate,
        so this loads a reaction SMILES through Ketcher's own API, the same way
        the page's "Example" button does (static/js/sketcher/integrated_ketcher.js). Ketcher then raises its change event,
        so the page reacts exactly as it does to a drawn reaction.

        Args:
            reaction_smiles: Reaction SMILES, e.g. "C.CC>>CP".
        """
        # the page hides its spinner once Ketcher is loaded and the reaction reloaded
        self.page.locator(".loading-circle").wait_for(
            state="hidden", timeout=KETCHER_LOAD_TIMEOUT_MS
        )
        # The page only starts listening for sketch changes a moment after the
        # spinner goes (static/js/sketcher/sketcher.js), with nothing on screen
        # to show it. So check the sketch reached the server, and draw it again
        # if it was set before the page was listening.
        for _ in range(DRAW_ATTEMPTS):
            try:
                with self.page.expect_request(
                    lambda request: request.url.endswith("/autoupdate_reaction_table"),
                    timeout=TABLE_REQUEST_TIMEOUT_MS,
                ):
                    self.page.evaluate(
                        "smiles => document.getElementById('ketcher-editor')"
                        ".contentWindow.ketcher.setMolecule(smiles)",
                        reaction_smiles,
                    )
                return
            except PlaywrightTimeoutError:
                continue
        raise AssertionError(
            f"Drawing {reaction_smiles} never asked the server for a reaction table"
        )

    def wait_for_reaction_table(self) -> None:
        """
        Wait for the server-rendered reaction table to appear.

        Building it takes a few seconds, longer than the usual assertion
        timeout, so wait for its first reactant before checking its contents.
        """
        self.reactant_name(1).wait_for(timeout=REACTION_TABLE_TIMEOUT_MS)

    @property
    def reaction_description(self) -> Locator:
        """The reaction description box at the top of the reaction table."""
        return self.page.get_by_placeholder("Please describe your reaction")

    def describe_reaction(self, description: str) -> None:
        """
        Enter the reaction's description and wait for the table to autosave.

        The table saves when an input changes, i.e. when it loses focus. Its
        contents are only kept for reopening the reaction once it has saved.

        Args:
            description: Text describing the reaction.
        """
        with self.page.expect_response(
            lambda response: response.url.endswith("/_autosave")
        ) as save:
            self.reaction_description.fill(description)
            self.reaction_description.blur()
        assert save.value.ok, f"Autosave failed with status {save.value.status}"

    def reactant_name(self, number: int) -> Locator:
        """
        Find a reactant's name in the reaction table.

        The table's inputs have no labels, so they are found by id.

        Args:
            number: The reactant's row number, starting at 1.

        Returns:
            Locator for the input holding the name.
        """
        return self.page.locator(f"#js-reactant{number}")

    def product_name(self, number: int) -> Locator:
        """
        Find a product's name in the reaction table.

        Args:
            number: The product's row number, starting at 1.

        Returns:
            Locator for the input holding the name.
        """
        return self.page.locator(f"#js-product{number}")

    @property
    def coshh_root(self) -> Locator:
        """The COSHH form section's body, where the React COSHH form is mounted."""
        return self.page.locator("#coshh-form-root")

    def open_coshh_section(self) -> None:
        """Expand the collapsed "COSHH Form" section."""
        self.page.get_by_role("button", name="COSHH Form").click()

    def create_coshh_form(self, template_name: str) -> None:
        """
        Create the reaction's COSHH form from a template. The section must be open.

        Args:
            template_name: Name of the template to base the form on.
        """
        choose_template(self.coshh_root, template_name)
        self.coshh_root.get_by_role("button", name="Create").click()

    @property
    def coshh_form(self) -> CoshhForm:
        """The reaction's COSHH form, once created."""
        return CoshhForm(self.page, self.coshh_root)
