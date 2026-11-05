"""Names of the seeded users and work structures that the end-to-end tests rely on."""

# from tests/populate_database.py; the main user is the Test-Workgroup PI
TEST_USERNAME = "test_username"
TEST_PASSWORD = "test_pw"
TEST_WORKGROUP = "Test-Workgroup"
TEST_WORKBOOK = "Test-Workbook"

# added by the e2e conftest to test workgroup role permissions
ROLES_WORKGROUP = "E2E-Roles-Workgroup"
SENIOR_RESEARCHER_USERNAME = "e2e_senior"
STANDARD_MEMBER_USERNAME = "e2e_member"
ROLE_USER_PASSWORD = "e2e_pw"
