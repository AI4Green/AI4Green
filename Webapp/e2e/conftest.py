"""
Session fixtures for the Playwright end-to-end suite.

The Flask app is built with the "test" config (so audit messages go to the
LoggingQueueProducer rather than Azure), seeded with the same data as the unit
tests, and served from a background thread on a free local port. Playwright's
``base_url`` points at that server, so tests navigate with relative paths.
"""

import threading
import uuid
from pathlib import Path
from typing import Callable, Iterator, List

import pytest
from e2e import populate_database as e2e_populate_database
from e2e.pages.home_page import HomePage
from e2e.pages.login_page import LoginPage
from e2e.seed_data import (
    ROLE_USER_PASSWORD,
    ROLES_WORKGROUP,
    SENIOR_RESEARCHER_USERNAME,
    STANDARD_MEMBER_USERNAME,
    TEST_PASSWORD,
    TEST_USERNAME,
)
from flask import Flask
from playwright.sync_api import Browser, BrowserContext, Page, expect
from sources import create_app
from sources.extensions import db
from tests import populate_database
from tests.utils import add_workgroup, main_test_person
from werkzeug.serving import BaseWSGIServer, make_server


def add_roles_workgroup() -> None:
    """
    Add a workgroup with one member in each role, for permission tests.

    The main test user is its principal investigator, plus a senior researcher
    and a standard member who are verified so they can log in.
    """
    senior_researcher = populate_database.add_test_user(
        username=SENIOR_RESEARCHER_USERNAME,
        email=f"{SENIOR_RESEARCHER_USERNAME}@test.com",
        fullname="E2E Senior Researcher",
        password=ROLE_USER_PASSWORD,
        verified=True,
    )
    standard_member = populate_database.add_test_user(
        username=STANDARD_MEMBER_USERNAME,
        email=f"{STANDARD_MEMBER_USERNAME}@test.com",
        fullname="E2E Standard Member",
        password=ROLE_USER_PASSWORD,
        verified=True,
    )
    add_workgroup(
        ROLES_WORKGROUP,
        principal_investigator=[main_test_person()],
        senior_researcher=[senior_researcher],
        standard_member=[standard_member],
    )


@pytest.fixture(scope="session")
def app() -> Iterator[Flask]:
    """
    Create the Flask app with a freshly created and seeded test database.

    Yields:
        The Flask application.
    """
    app = create_app("test")

    with app.app_context():
        db.drop_all()
        db.create_all()
        # adding users records a location, which is read from the request json
        with app.test_request_context(
            json={
                "IP_address": "104.28.100.0",
                "country": "United Kingdom",
                "city": "Nottingham",
            }
        ):
            populate_database.insert_test_data()
            add_roles_workgroup()
            e2e_populate_database.add_coshh_data()
            e2e_populate_database.add_inventory_data()

    yield app

    with app.app_context():
        db.session.remove()
        db.drop_all()


@pytest.fixture(scope="session")
def live_server(app: Flask) -> Iterator[str]:
    """
    Serve the app from a background thread on a free port.

    Args:
        app: The seeded Flask application.

    Yields:
        The root URL of the running server.
    """
    server: BaseWSGIServer = make_server("127.0.0.1", 0, app, threaded=True)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    yield f"http://127.0.0.1:{server.server_port}"

    server.shutdown()
    thread.join()


@pytest.fixture(scope="session")
def base_url(live_server: str) -> str:
    """
    Override pytest-base-url so Playwright resolves relative URLs against the live server.

    Args:
        live_server: The root URL of the running server.

    Returns:
        The base URL used by every browser context.
    """
    return live_server


def track_page_errors(page: Page) -> List[str]:
    """
    Record uncaught JavaScript exceptions raised on a page.

    Args:
        page: The Playwright page to listen to.

    Returns:
        A list that fills with error messages and stack traces as they happen.
    """
    errors: List[str] = []
    page.on("pageerror", lambda error: errors.append(error.stack or error.message))
    return errors


@pytest.fixture
def page(page: Page) -> Iterator[Page]:
    """
    Wrap pytest-playwright's ``page`` so any uncaught JavaScript error fails the test.

    Args:
        page: The page provided by pytest-playwright.

    Yields:
        The same page, with error tracking attached.
    """
    errors = track_page_errors(page)
    yield page
    assert not errors, f"Uncaught JavaScript errors on the page: {errors}"


@pytest.fixture(scope="session")
def logged_in_storage_state(
    browser: Browser, base_url: str, tmp_path_factory: pytest.TempPathFactory
) -> Path:
    """
    Log in once through the UI and save the session cookie for reuse.

    Args:
        browser: The session-wide Playwright browser.
        base_url: The live server URL.
        tmp_path_factory: Pytest factory for temporary directories.

    Returns:
        Path to the saved storage state JSON.
    """
    context = browser.new_context(base_url=base_url)
    login_page = LoginPage(context.new_page())
    login_page.goto()
    login_page.login(TEST_USERNAME, TEST_PASSWORD)
    expect(HomePage(login_page.page).welcome_heading).to_be_visible()

    storage_state_path = tmp_path_factory.mktemp("auth") / "storage_state.json"
    context.storage_state(path=storage_state_path)
    context.close()
    return storage_state_path


@pytest.fixture
def logged_in_page(new_context, logged_in_storage_state: Path) -> Iterator[Page]:
    """
    Provide a page that is already logged in as the seeded test user.

    Uses pytest-playwright's ``new_context`` so the CLI tracing and screenshot
    options still apply.

    Args:
        new_context: pytest-playwright factory for browser contexts.
        logged_in_storage_state: Saved cookies for the logged-in user.

    Yields:
        A logged-in page with error tracking attached.
    """
    context: BrowserContext = new_context(storage_state=logged_in_storage_state)
    logged_in_page = context.new_page()
    errors = track_page_errors(logged_in_page)
    yield logged_in_page
    assert not errors, f"Uncaught JavaScript errors on the page: {errors}"


@pytest.fixture
def log_in_as(new_context) -> Iterator[Callable[[str, str], Page]]:
    """
    Provide a function that logs in as any seeded user in a fresh browser context.

    Use this for users other than the main test user, e.g. to check role permissions.

    Args:
        new_context: pytest-playwright factory for browser contexts.

    Yields:
        A function taking a username and password that returns a logged-in page.
    """
    errors_per_page: List[List[str]] = []

    def _log_in_as(username: str, password: str) -> Page:
        context: BrowserContext = new_context()
        page = context.new_page()
        errors_per_page.append(track_page_errors(page))
        login_page = LoginPage(page)
        login_page.goto()
        login_page.login(username, password)
        expect(HomePage(page).welcome_heading).to_contain_text(username)
        return page

    yield _log_in_as
    errors = [error for page_errors in errors_per_page for error in page_errors]
    assert not errors, f"Uncaught JavaScript errors on the page: {errors}"


@pytest.fixture
def inventory_container(app: Flask) -> str:
    """
    Add a fresh container of the stocked compound for one test to act on.

    Args:
        app: The seeded Flask application.

    Returns:
        The container's lot number, which is unique and labels its row actions.
    """
    lot_number = f"E2E-{uuid.uuid4().hex[:8].upper()}"
    with app.app_context():
        e2e_populate_database.add_inventory_container(lot_number)
    return lot_number
