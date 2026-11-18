"""Names of the seeded users and work structures that the end-to-end tests rely on."""

# from tests/populate_database.py; the main user is the Test-Workgroup PI
TEST_USERNAME = "test_username"
TEST_PASSWORD = "test_pw"
TEST_WORKGROUP = "Test-Workgroup"
TEST_WORKBOOK = "Test-Workbook"
# a reaction already in Test-Workbook, so new reactions can't reuse its name
SEEDED_REACTION_NAME = "a test reaction"

# added by the e2e conftest to test workgroup role permissions
ROLES_WORKGROUP = "E2E-Roles-Workgroup"
SENIOR_RESEARCHER_USERNAME = "e2e_senior"
STANDARD_MEMBER_USERNAME = "e2e_member"
ROLE_USER_PASSWORD = "e2e_pw"

# added by e2e/populate_database.py for the COSHH tests
HS_ADMIN_USERNAME = "e2e_hs_admin"
COSHH_TEMPLATE = "E2E COSHH Template"
COSHH_SECTION = "E2E Hazard Assessment"
COSHH_TEXT_FIELD = "Control measures"

# added by e2e/populate_database.py for the inventory tests, in ROLES_WORKGROUP
INVENTORY_WORKGROUP = ROLES_WORKGROUP
INVENTORY_LOCATION = "E2E Shelf A"
INVENTORY_OTHER_LOCATION = "E2E Shelf B"
# compounds from tests/populate_database.py: one already stocked, one not
STOCKED_COMPOUND_NAME = "Testoic Acid"
STOCKED_COMPOUND_CAS = "123-45-5"
UNSTOCKED_COMPOUND_NAME = "Testide"
UNSTOCKED_COMPOUND_CAS = "123-45-8"
