#
# Shroud tests
#
# (from parent directory) python -m unittest test
#
from __future__ import absolute_import

import unittest

from . import test_parse_decl
from . import test_shroud
from . import test_util


test_cases = (
    test_parse_decl.CheckDeclCase,
    test_util.UtilCase,
    test_util.OptionCase,
    test_shroud.MainCase,
)


def load_tests(loader, tests, pattern):
    suite = unittest.TestSuite()
    for test_class in test_cases:
        tests = loader.loadTestsFromTestCase(test_class)
        suite.addTests(tests)
    return suite

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(load_tests())
