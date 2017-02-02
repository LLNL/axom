from __future__ import print_function
#from __future__ import unicode_literals

from shroud import main

import os
import sys
import sysconfig
import unittest

import do_test

def distutils_dir_name(dname):
    """Returns the name of a distutils build directory"""
    f = "{dirname}.{platform}-{version[0]}.{version[1]}"
    return f.format(dirname=dname,
                    platform=sysconfig.get_platform(),
                    version=sys.version_info)


class MainCase(unittest.TestCase):
    def setUp(self):

        # python -m unittest tests

#        self.cwd = os.getcwd()
        self.cwd = os.path.join(os.getcwd(), 'tests')

        self.testdir = os.path.abspath(
            os.path.join(
                self.cwd, '..', 'build', distutils_dir_name('temp')))
        self.tester = do_test.Tester()
        self.assertTrue(
            self.tester.set_environment(self.cwd, self.testdir))

    def run_shroud(self, input):
        tester = self.tester
        tester.open_log(input + '.log')
        tester.set_test(input)
        status = tester.do_module()
        tester.close_log()
        self.assertTrue(status)

    def test_example(self):
        self.run_shroud('example')

    def test_include(self):
        self.run_shroud('include')

    def test_names(self):
        self.run_shroud('names')

    def test_strings(self):
        self.run_shroud('strings')

    def test_tutorial(self):
        self.run_shroud('tutorial')


if __name__ == '__main__':
    unittest.main()
