from __future__ import print_function

from shroud import main

import argparse
import io
import sys
import unittest


class MainCase(unittest.TestCase):
    def setUp(self):
        # Default command line arguments
        self.args = argparse.Namespace(
            outdir='',
            outdir_c_fortran='',
            outdir_python='',
            outidr_lua='',
            logdir='',
            cfiles='',
            ffiles='',
            path=[],
            filename=[]
        )

        # redirect stdout and stderr
        self.stdout = io.StringIO()
        self.saved_stdout = sys.stdout
        sys.stdout = self.stdout

        self.stderr = io.StringIO()
        self.saved_stderr = sys.stderr
        sys.stdout = self.stderr

    def tearDown(self):
        self.stdout.close()
        sys.stdout = self.saved_stdout

        self.stderr.close()
        sys.stderr = self.saved_stderr

    def test_no_args(self):
        """Run with no arguments.
        """
        args = self.args
        with self.assertRaises(SystemExit):
            main.main_with_args(args)


if __name__ == '__main__':
    unittest.main()
