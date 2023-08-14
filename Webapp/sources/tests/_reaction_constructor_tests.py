import test_accessibility
import test_add_reagent_new_cas
import test_complete_reaction_locked
import test_compound_data_error_report
import test_demo
import test_novel_compounds
import test_novel_solvents
import test_pdf_summary
import test_primary_product_change
import test_reaction_addenda
import test_reaction_attachments
import test_reactant_tickbox_mass
import test_reaction_reload
import test_reaction_table
import test_reaction_table_button_hovers
import test_reaction_table_missing_data
import test_reaction_table_numbering
import test_reaction_tutorial
import test_reagent_list
import test_reagent_selecting
import test_remove_any_order
import test_resubmit_dataloss_warning
import test_save_reaction
import test_solvent_guide
import test_solvent_selection
import test_summary_rendered
import test_summary_table
import test_unknown_reactant_product
import test_wb_reactions_from_all


import unittest


def reaction_suite_test():
    loader = unittest.TestLoader()
    suite_test = unittest.TestSuite()
    suite_test.addTests(loader.loadTestsFromModule(test_accessibility))
    suite_test.addTests(loader.loadTestsFromModule(test_add_reagent_new_cas))
    suite_test.addTests(loader.loadTestsFromModule(test_complete_reaction_locked))
    suite_test.addTests(loader.loadTestsFromModule(test_compound_data_error_report))
    suite_test.addTests(loader.loadTestsFromModule(test_demo))
    suite_test.addTests(loader.loadTestsFromModule(test_novel_compounds))
    suite_test.addTests(loader.loadTestsFromModule(test_novel_solvents))
    suite_test.addTests(loader.loadTestsFromModule(test_pdf_summary))
    suite_test.addTests(loader.loadTestsFromModule(test_primary_product_change))
    suite_test.addTests(loader.loadTestsFromModule(test_reactant_tickbox_mass))
    suite_test.addTests(loader.loadTestsFromModule(test_reaction_addenda))
    suite_test.addTests(loader.loadTestsFromModule(test_reaction_attachments))
    suite_test.addTests(loader.loadTestsFromModule(test_reaction_reload))
    suite_test.addTests(loader.loadTestsFromModule(test_reaction_table))
    suite_test.addTests(loader.loadTestsFromModule(test_reaction_table_button_hovers))
    suite_test.addTests(loader.loadTestsFromModule(test_reaction_table_missing_data))
    suite_test.addTests(loader.loadTestsFromModule(test_reaction_table_numbering))
    suite_test.addTests(loader.loadTestsFromModule(test_reaction_tutorial))
    suite_test.addTests(loader.loadTestsFromModule(test_reagent_list))
    suite_test.addTests(loader.loadTestsFromModule(test_reagent_selecting))
    suite_test.addTests(loader.loadTestsFromModule(test_resubmit_dataloss_warning))
    suite_test.addTests(loader.loadTestsFromModule(test_remove_any_order))
    suite_test.addTests(loader.loadTestsFromModule(test_save_reaction))
    suite_test.addTests(loader.loadTestsFromModule(test_solvent_guide))
    suite_test.addTests(loader.loadTestsFromModule(test_solvent_selection))
    suite_test.addTests(loader.loadTestsFromModule(test_summary_rendered))
    suite_test.addTests(loader.loadTestsFromModule(test_summary_table))
    suite_test.addTests(loader.loadTestsFromModule(test_unknown_reactant_product))
    suite_test.addTests(loader.loadTestsFromModule(test_wb_reactions_from_all))
    return suite_test


if __name__ == '__main__':

    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(reaction_suite_test())
