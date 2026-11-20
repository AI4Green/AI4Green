"""Page object for the workup and purification cards on a reaction's page."""

import re
from contextlib import contextmanager
from typing import Iterator

from e2e.pages.base_page import BasePage
from playwright.sync_api import Locator, Page

# the step card saves 800 ms after the last keystroke; allow time for the request
STEP_SAVE_TIMEOUT_MS = 10_000


class WorkupSection(BasePage):
    """
    One of the React workup cards in the reaction page's "Reaction Notes"
    section: "Reaction Workup" or "Reaction Purification".
    """

    def __init__(self, page: Page, category: str = "WORKUP") -> None:
        """
        Args:
            page: The Playwright page to drive.
            category: "WORKUP" or "PURIFICATION".
        """
        super().__init__(page)
        list_id = (
            "#reaction-workup-list"
            if category == "WORKUP"
            else "#reaction-purification-list"
        )
        self.root = page.locator(list_id)

    def open(self) -> None:
        """Expand the collapsed "Reaction Notes" section, which holds the cards."""
        self.page.get_by_role("button", name="Reaction Notes").click()
        self.root.get_by_role("button", name="Add step").wait_for()

    def add_step(self, step_label: str) -> Locator:
        """
        Add a step from the "Add step" menu.

        Args:
            step_label: The step type's label in the menu, e.g. "Wash".

        Returns:
            The new step's card.
        """
        count = self.steps.count()
        self.root.get_by_role("button", name="Add step").click()
        # Chakra renders the menu in a portal at the end of <body>
        self.page.get_by_role("menuitem", name=step_label, exact=True).click()
        return self.step(count + 1, step_label)

    def step(self, position: int, step_label: str) -> Locator:
        """
        Find a step's card.

        Args:
            position: The step's position, starting at 1.
            step_label: The step type's label, e.g. "Wash".

        Returns:
            Locator for the card.
        """
        return self.root.get_by_role("region", name=f"Step {position}: {step_label}")

    @property
    def steps(self) -> Locator:
        """Every step card in this card, in order."""
        return self.root.get_by_role("region", name=re.compile(r"^Step \d+: "))

    def expand_step(self, step: Locator) -> None:
        """
        Expand a step's card to show its inputs. Cards start collapsed, showing
        only a one-line summary.

        Args:
            step: The step's card.
        """
        step.get_by_role("button", name="Expand step").click()

    @contextmanager
    def _saving_step(self) -> Iterator[None]:
        """Wait for the step autosave that the wrapped action triggers, and check it worked."""
        with self.page.expect_response(
            lambda response: re.search(r"/api/workups/.+/steps/\d+$", response.url)
            is not None
            and response.request.method == "PUT",
            timeout=STEP_SAVE_TIMEOUT_MS,
        ) as save:
            yield
        assert save.value.ok, f"Step save failed with status {save.value.status}"

    def fill_and_wait_for_save(self, step: Locator, label: str, value: str) -> None:
        """
        Type into one of a step's fields and wait for the step to autosave.

        Args:
            step: The step's card, expanded.
            label: The field's label, e.g. "Organic phase".
            value: Text to enter.
        """
        with self._saving_step():
            step.get_by_label(label, exact=True).fill(value)

    def select_and_wait_for_save(self, step: Locator, label: str, option: str) -> None:
        """
        Choose an option in one of a step's dropdowns and wait for the step to
        autosave.

        Args:
            step: The step's card, expanded.
            label: The dropdown's label, e.g. "Method".
            option: The option to choose, e.g. "Rotary evaporator".
        """
        with self._saving_step():
            step.get_by_label(label, exact=True).select_option(option)

    def choose_solvent_and_wait_for_save(
        self, step: Locator, label: str, search: str, solvent: str
    ) -> None:
        """
        Pick a solvent from a material's coloured CHEM21 solvent list and wait for
        the step to autosave.

        Args:
            step: The step's card, expanded.
            label: The material's label, e.g. "Organic phase".
            search: Text to type to narrow the list.
            solvent: The solvent to pick, e.g. "Testanol".
        """
        step.get_by_label(label, exact=True).fill(search)
        with self._saving_step():
            step.get_by_role(
                "option", name=re.compile(f"^{re.escape(solvent)}")
            ).click()

    def save_as_template(self, name: str) -> None:
        """
        Save the card's steps as a workbook template.

        Args:
            name: The template's name.
        """
        self.root.get_by_role("button", name="Save as template").click()
        dialog = self.page.get_by_role("dialog", name="Save workup as template")
        dialog.get_by_label("Template name").fill(name)
        dialog.get_by_role("button", name="Save template").click()

    def start_from_template(self, name: str) -> None:
        """
        Add the steps of a workbook template.

        Args:
            name: The template's name.
        """
        self.root.get_by_role("button", name="Start from template").click()
        dialog = self.page.get_by_role("dialog", name="workup templates")
        dialog.get_by_role("listitem", name=name).get_by_role(
            "button", name="Use"
        ).click()
