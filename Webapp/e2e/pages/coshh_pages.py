"""
Page objects for COSHH forms in the React app.

H&S admins design COSHH templates on COSHH Form Management. Users fill in COSHH
forms made from those templates, either standalone from My COSHH Forms or
embedded in a reaction's page. The form itself is the same component in both
places, so ``CoshhForm`` is shared.
"""

from e2e.pages.base_page import BasePage
from playwright.sync_api import Locator, Page

# TextAreaField (client-app/src/components/core/forms/textarea-field.jsx)
# debounces by 150 ms; wait a little longer than that
TEXTAREA_DEBOUNCE_WAIT_MS = 400


def choose_template(container: Locator, template_name: str) -> None:
    """
    Pick a template in the "COSHH Template" dropdown, a react-select combobox.

    Args:
        container: Element holding the dropdown, e.g. the create modal.
        template_name: Name of the template to pick.
    """
    combobox = container.get_by_role("combobox")
    combobox.click()
    combobox.fill(template_name)
    combobox.press("Enter")


class CoshhForm:
    """A COSHH form being filled in: one block per template section, which autosaves."""

    def __init__(self, page: Page, container: Locator) -> None:
        """
        Args:
            page: The Playwright page to drive.
            container: Element the form is rendered in.
        """
        self.page = page
        self.container = container

    def section(self, section_name: str) -> Locator:
        """
        Find a section's title.

        Args:
            section_name: The section's name in the template.

        Returns:
            Locator for the title text.
        """
        return self.container.get_by_text(section_name, exact=True)

    def text_field(self, field_name: str) -> Locator:
        """
        Find a text field by its label.

        Args:
            field_name: The field's name in the template.

        Returns:
            Locator for the input.
        """
        return self.container.get_by_label(field_name)

    def fill_and_wait_for_autosave(self, field_name: str, value: str) -> None:
        """
        Type into a text field and wait for the autosave it triggers to succeed.

        The form saves itself about a second after the last change, with no
        message on screen, so this waits for the save request instead.

        Args:
            field_name: The field's name in the template.
            value: Text to enter.
        """
        with self.page.expect_response(
            lambda response: "/api/template-instances/" in response.url
            and response.request.method == "PUT"
        ) as save:
            self.text_field(field_name).fill(value)
        assert save.value.ok, f"Autosave failed with status {save.value.status}"


class CoshhFormManagementPage(BasePage):
    """Template list at /spa/coshh-form-management, for designing COSHH templates."""

    PATH = "/spa/coshh-form-management"

    def goto(self) -> None:
        """Open COSHH Form Management."""
        self.page.goto(self.PATH)

    @property
    def new_button(self) -> Locator:
        """The "New" button, shown only to H&S admins."""
        return self.page.get_by_role("button", name="New", exact=True)

    @property
    def create_dialog(self) -> Locator:
        """The "Create Project Type" modal for a new template, while open."""
        return self.page.get_by_role("dialog", name="Create Project Type")

    def create_template(self, name: str, description: str) -> None:
        """
        Create an empty template with "New".

        Args:
            name: Template name, unique among templates.
            description: Template description, required.
        """
        self.new_button.click()
        dialog = self.create_dialog
        dialog.get_by_label("Project type name").fill(name)
        dialog.get_by_label("Project type description").fill(description)
        # TextAreaField only passes its value to the form after a debounce, with
        # nothing on screen to wait for, so give it time before submitting
        self.page.wait_for_timeout(TEXTAREA_DEBOUNCE_WAIT_MS)
        dialog.get_by_role("button", name="Create").click()

    def template_row(self, name: str) -> Locator:
        """
        Find a template's row in the list.

        Args:
            name: The template's name.

        Returns:
            Locator for the row.
        """
        return self.page.get_by_role("row").filter(
            has=self.page.get_by_role("cell", name=name, exact=True)
        )


class MyCoshhFormsPage(BasePage):
    """The current user's COSHH forms at /spa/my-coshh-forms."""

    PATH = "/spa/my-coshh-forms"

    def goto(self) -> None:
        """Open My COSHH Forms."""
        self.page.goto(self.PATH)

    def create_form(self, template_name: str) -> None:
        """
        Start a standalone COSHH form from a template with "New".

        On success the app opens the new form.

        Args:
            template_name: Name of the template to base the form on.
        """
        self.page.get_by_role("button", name="New", exact=True).click()
        dialog = self.page.get_by_role("dialog", name="Create COSHH Form")
        choose_template(dialog, template_name)
        dialog.get_by_role("button", name="Create").click()

    @property
    def form(self) -> CoshhForm:
        """The open COSHH form."""
        # the React app's root element
        return CoshhForm(self.page, self.page.locator("#root"))

    def form_rows(self, template_name: str) -> Locator:
        """
        Find the listed forms made from a template.

        Args:
            template_name: Name of the template.

        Returns:
            Locator for the matching rows.
        """
        return self.page.get_by_role("row").filter(
            has=self.page.get_by_role("cell", name=template_name, exact=True)
        )
