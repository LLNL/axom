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
    # used from 'python -m unittest tests'
    suite = unittest.TestSuite()
    for test_class in test_cases:
        tests = loader.loadTestsFromTestCase(test_class)
        suite.addTests(tests)
    return suite


def load_tests2():
    # used from 'setup.py test'
    loader = unittest.TestLoader()
    return load_tests(loader, None, None)

