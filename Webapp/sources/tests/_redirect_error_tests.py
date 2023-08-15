import test_reagents
import test_solvent_cas
import test_solvents

import unittest


def redirect_error_suite_test():
	"""All these tests fail on redirect issues when ran alongside other tests
	Removed from reaction constructor tests and __init__
	"""
	loader = unittest.TestLoader()
	suite_test = unittest.TestSuite()
	suite_test.addTests(loader.loadTestsFromModule(test_reagents))
	suite_test.addTests(loader.loadTestsFromModule(test_solvent_cas))
	suite_test.addTests(loader.loadTestsFromModule(test_solvents))
	return suite_test

if __name__ == '__main__':
	runner = unittest.TextTestRunner(verbosity=2)
	runner.run(redirect_error_suite_test())
