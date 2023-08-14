import test_accessibility
import test_about_page
import test_add_reagent_new_cas
import test_admin_dashboard
import test_auth
import test_buttons_rendered_linked
import test_complete_reaction_locked
import test_compound_data_error_report
import test_cookies
# import test_input_reagent
import test_demo
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
import test_novel_compounds
import test_novel_solvents
import test_pdf_summary
import test_pi_sr_request
import test_primary_product_change
import test_reactant_tickbox_mass
import test_reaction_addenda
import test_reaction_attachments
import test_reaction_list
import test_reaction_reload
import test_reaction_table
import test_reaction_table_button_hovers
import test_reaction_table_missing_data
import test_reaction_table_numbering
import test_reaction_tutorial
import test_reagent_list
import test_reagent_selecting
import test_reagents
import test_remove_any_order
import test_reset_password
import test_resubmit_dataloss_warning
import test_return_to_login_page
import test_save_reaction
import test_search
import test_solvent_cas
import test_solvent_guide
import test_solvent_selection
import test_solvents
import test_summary_rendered
import test_summary_table
import test_unknown_reactant_product
import test_user_database
import test_wb_reactions_from_all
import test_WG_summary
import test_workbook_creation
import test_workgroup
import test_workgroup_creation

import unittest

print('main tests')
# Edit run configuration to run with environmental variables: TESTING=1

def suite_test():
    loader = unittest.TestLoader()
    suite_test = unittest.TestSuite()
    suite_test.addTests(loader.loadTestsFromModule(test_accessibility))
    suite_test.addTests(loader.loadTestsFromModule(test_about_page))
    suite_test.addTests(loader.loadTestsFromModule(test_add_reagent_new_cas))
    suite_test.addTests(loader.loadTestsFromModule(test_admin_dashboard))
    suite_test.addTests(loader.loadTestsFromModule(test_auth))
    suite_test.addTests(loader.loadTestsFromModule(test_buttons_rendered_linked))
    suite_test.addTests(loader.loadTestsFromModule(test_complete_reaction_locked))
    suite_test.addTests(loader.loadTestsFromModule(test_compound_data_error_report))
    suite_test.addTests(loader.loadTestsFromModule(test_cookies))
    # suite_test.addTests(loader.loadTestsFromModule(test_input_reagent))
    suite_test.addTests(loader.loadTestsFromModule(test_demo))
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
    suite_test.addTests(loader.loadTestsFromModule(test_novel_compounds))
    suite_test.addTests(loader.loadTestsFromModule(test_novel_solvents))
    suite_test.addTests(loader.loadTestsFromModule(test_pdf_summary))
    suite_test.addTests(loader.loadTestsFromModule(test_pi_sr_request))
    suite_test.addTests(loader.loadTestsFromModule(test_primary_product_change))
    suite_test.addTests(loader.loadTestsFromModule(test_reactant_tickbox_mass))
    suite_test.addTests(loader.loadTestsFromModule(test_reaction_addenda))
    suite_test.addTests(loader.loadTestsFromModule(test_reaction_attachments))
    suite_test.addTests(loader.loadTestsFromModule(test_reaction_list))
    suite_test.addTests(loader.loadTestsFromModule(test_reaction_reload))
    suite_test.addTests(loader.loadTestsFromModule(test_reaction_table))
    suite_test.addTests(loader.loadTestsFromModule(test_reaction_table_button_hovers))
    suite_test.addTests(loader.loadTestsFromModule(test_reaction_table_missing_data))
    suite_test.addTests(loader.loadTestsFromModule(test_reaction_table_numbering))
    suite_test.addTests(loader.loadTestsFromModule(test_reaction_tutorial))
    suite_test.addTests(loader.loadTestsFromModule(test_reagent_list))
    suite_test.addTests(loader.loadTestsFromModule(test_reagent_selecting))
    # suite_test.addTests(loader.loadTestsFromModule(test_reagents))
    suite_test.addTests(loader.loadTestsFromModule(test_reset_password))
    suite_test.addTests(loader.loadTestsFromModule(test_resubmit_dataloss_warning))
    # suite_test.addTests(loader.loadTestsFromModule(test_remove_from_WB_WG))
    suite_test.addTests(loader.loadTestsFromModule(test_remove_any_order))
    suite_test.addTests(loader.loadTestsFromModule(test_return_to_login_page))
    suite_test.addTests(loader.loadTestsFromModule(test_save_reaction))
    suite_test.addTests(loader.loadTestsFromModule(test_search))
    # suite_test.addTests(loader.loadTestsFromModule(test_solvent_cas))
    suite_test.addTests(loader.loadTestsFromModule(test_solvent_guide))
    suite_test.addTests(loader.loadTestsFromModule(test_solvent_selection))
    # suite_test.addTests(loader.loadTestsFromModule(test_solvents))
    suite_test.addTests(loader.loadTestsFromModule(test_summary_rendered))
    suite_test.addTests(loader.loadTestsFromModule(test_summary_table))
    suite_test.addTests(loader.loadTestsFromModule(test_unknown_reactant_product))
    suite_test.addTests(loader.loadTestsFromModule(test_user_database))
    suite_test.addTests(loader.loadTestsFromModule(test_wb_reactions_from_all))
    suite_test.addTests(loader.loadTestsFromModule(test_WG_summary))
    suite_test.addTests(loader.loadTestsFromModule(test_workbook_creation))
    suite_test.addTests(loader.loadTestsFromModule(test_workgroup))
    suite_test.addTests(loader.loadTestsFromModule(test_workgroup_creation))
    return suite_test


if __name__ == '__main__':

    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite_test())
