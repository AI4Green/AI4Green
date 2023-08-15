import test_accessibility
import test_admin_dashboard
import test_auth
import test_buttons_rendered_linked
import test_cookies
import test_help_pages
import test_join_workbook
import test_join_workgroup
import test_landing_page
import test_manage_account
import test_manage_workbook
import test_manage_workgroup
import test_news_feed
import test_notification_email
import test_notification_read
import test_search
import test_pi_sr_request
import test_reaction_list
import test_reset_password
import test_return_to_login_page
import test_user_database
import test_wb_reactions_from_all
import test_WG_summary
import test_workbook_creation
import test_workgroup
import test_workgroup_creation

import unittest


def non_reaction_suite_test():
    loader = unittest.TestLoader()
    suite_test = unittest.TestSuite()
    suite_test.addTests(loader.loadTestsFromModule(test_accessibility))
    suite_test.addTests(loader.loadTestsFromModule(test_admin_dashboard))
    suite_test.addTests(loader.loadTestsFromModule(test_auth))
    suite_test.addTests(loader.loadTestsFromModule(test_buttons_rendered_linked))
    suite_test.addTests(loader.loadTestsFromModule(test_cookies))
    suite_test.addTests(loader.loadTestsFromModule(test_help_pages))
    suite_test.addTests(loader.loadTestsFromModule(test_join_workbook))
    suite_test.addTests(loader.loadTestsFromModule(test_join_workgroup))
    suite_test.addTests(loader.loadTestsFromModule(test_landing_page))
    suite_test.addTests(loader.loadTestsFromModule(test_manage_account))
    suite_test.addTests(loader.loadTestsFromModule(test_manage_workbook))
    suite_test.addTests(loader.loadTestsFromModule(test_manage_workgroup))
    suite_test.addTests(loader.loadTestsFromModule(test_news_feed))
    suite_test.addTests(loader.loadTestsFromModule(test_notification_email))
    suite_test.addTests(loader.loadTestsFromModule(test_notification_read))
    suite_test.addTests(loader.loadTestsFromModule(test_search))
    suite_test.addTests(loader.loadTestsFromModule(test_pi_sr_request))
    suite_test.addTests(loader.loadTestsFromModule(test_reaction_list))
    suite_test.addTests(loader.loadTestsFromModule(test_reset_password))
    suite_test.addTests(loader.loadTestsFromModule(test_return_to_login_page))
    suite_test.addTests(loader.loadTestsFromModule(test_user_database))
    suite_test.addTests(loader.loadTestsFromModule(test_wb_reactions_from_all))
    suite_test.addTests(loader.loadTestsFromModule(test_WG_summary))
    suite_test.addTests(loader.loadTestsFromModule(test_workbook_creation))
    suite_test.addTests(loader.loadTestsFromModule(test_workgroup))
    suite_test.addTests(loader.loadTestsFromModule(test_workgroup_creation))
    return suite_test


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(non_reaction_suite_test())
