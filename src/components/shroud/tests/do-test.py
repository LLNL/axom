#!/usr/bin/env python
#
# Run tests for shroud
# run input file, compare results to reference files
#
# Directory structure
#  src-dir/name.yaml
#  src-dir/ref/name/      reference results
#
#  bin-dir/test/name/     generated files

# logging.debug('This is a debug message')
# logging.info('This is an info message')
# logging.warning('This is a warning message')
# logging.error('This is an error message')
# logging.critical('This is a critical error message')

from __future__ import print_function

import argparse
import errno
import filecmp
import logging
import os
import subprocess
import sys

#from io import StringIO
from io import BytesIO as StringIO

import shroud.main

#subprocess.call("ls -l", shell=True)

#proc = subprocess.Popen(['tail', '-500', 'mylogfile.log'], stdout=subprocess.PIPE)
#for line in proc.stdout.readlines():
#    print line.rstrip()


class Tester:
    def __init__(self):
        self.test_input_dir = ''
        self.test_output_dir = ''

        self.code_path = ''

        self.testyaml = ''   # input file
        self.ref_dir = ''    # reference directory
        self.result_dir = ''

    def open_log(self, logname):
        logging.basicConfig(filename=os.path.join(
            self.test_output_dir, logname),
                            filemode='w',
                            level=logging.DEBUG,
                        )

    def close_log(self):
        logging.shutdown()

    def set_environment(self, input, output, executable=None):
        """Set environment for all tests.
        """
        self.test_input_dir = input
        self.test_output_dir = output

        status = True
        if not os.path.isdir(input):
            status = False
            print('Missing source directory:', input)
        if executable:
            if not os.path.isdir(executable):
                status = False
                print('Missing executable directory:', executable)
            self.code_path = os.path.join(executable, 'shroud')
        makedirs(output)
        return status

    def set_test(self, name, replace_ref=False):
        """Setup for a single test.     
        """
        self.testname = name
        logging.info('--------------------------------------------------')
        logging.info('Testing ' + name)

        self.testyaml = os.path.join(self.test_input_dir, name + '.yaml')
        logging.info('Input file: ' + self.testyaml)
        if not os.path.isfile(self.testyaml):
            logging.error('Input file does not exist')
            return False

        self.ref_dir = os.path.join(self.test_input_dir, name)
        logging.info('Reference directory: ' + self.ref_dir)

        if replace_ref:
            # replacing reference, just create directly in ref directory
            self.result_dir = self.ref_dir
        else:
            self.result_dir = os.path.join(self.test_output_dir, name)
            logging.info('Result directory: ' + self.result_dir)
            makedirs(self.result_dir)
            clear_files(self.result_dir)

        return True

    def push_stdout(self):
        # redirect stdout and stderr
        self.stdout = StringIO()
        self.saved_stdout = sys.stdout
        sys.stdout = self.stdout

        self.stderr = StringIO()
        self.saved_stderr = sys.stderr
        sys.stderr = self.stderr

    def pop_stdout(self):
        self.stdout_lines = self.stdout.getvalue()
        self.stdout.close()
        sys.stdout = self.saved_stdout

        self.stderr_lines = self.stderr.getvalue()
        self.stderr.close()
        sys.stderr = self.saved_stderr

    def do_module(self):
        """Run Shroud via a method."""
        args = argparse.Namespace(
            outdir=self.result_dir,
            outdir_c_fortran='',
            outdir_python='',
            outdir_lua='',
            logdir=self.result_dir,
            cfiles='',
            ffiles='',
            path=[self.test_input_dir],
            filename=[self.testyaml]
        )
        logging.info('Arguments: ' + str(args))

        status = True
        self.push_stdout()
        try:
            shroud.main.main_with_args(args)
        except:
            logging.error('Shroud failed')
            status = False
        self.pop_stdout()

        # write output to a file
        output_file = os.path.join(self.result_dir, 'output')
        fp = open(output_file, 'w')
        fp.write(self.stdout_lines)
        fp.close()

        if status:
            status = self.do_compare()

        return status

    def do_test(self):
        """ Run test, return True/False for pass/fail.
        Files must compare, with no extra or missing files.
        """
        logging.info('Code to test: ' + self.code_path)

        cmd = [
            self.code_path,
            '--path', self.test_input_dir,
            '--logdir', self.result_dir,
            '--outdir', self.result_dir,
            self.testyaml,
            ]
        logging.debug(' '.join(cmd))

        try:
            output = subprocess.check_output(
                cmd,
                stderr=subprocess.STDOUT,
                universal_newlines=True)
        except subprocess.CalledProcessError as exc:
            logging.error('Exit status: %d' % exc.returncode)
            logging.error(exc.output)
            return False

        output_file = os.path.join(self.result_dir, 'output')
        fp = open(output_file, 'w')
        fp.write(output)
        fp.close()

        return True

    def do_compare(self):
        status = True  # assume it passes

        cmp = filecmp.dircmp(self.ref_dir, self.result_dir)
        if not os.path.exists(self.ref_dir):
            logging.info('Reference directory does not exist: ' + self.ref_dir)
            return False

        match, mismatch, errors = filecmp.cmpfiles(self.ref_dir, self.result_dir, cmp.common)
        for file in cmp.common:
            logging.info('Compare: ' + file)
        if mismatch:
            status = False
            for file in mismatch:
                logging.warn('Does not compare: '+ file)
        if errors:
            status = False
            for file in errors:
                logging.warn('Unable to compare: ' + file)

        if cmp.left_only:
            status = False
            for file in cmp.left_only:
                logging.warn('Only in reference: ' + file)
        if cmp.right_only:
            status = False
            for file in cmp.right_only:
                logging.warn('Only in result: ' + file)

        if status:
            logging.info('Test {} pass'.format(self.testname))
        else:
            logging.info('Test {} fail'.format(self.testname))
        return status


def makedirs(path):
    """ Make sure directory exists.
    """
    try:
        # Python >=3.2
        os.makedirs(path, exist_ok=True)
    except TypeError:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST or not os.path.isdir(path):
                raise


def clear_files(path):
    """Remove all files in a directory.
    """
    for file in os.listdir(path):
        full_path = os.path.join(path, file)
        try:
            if os.path.isfile(full_path):
                os.unlink(full_path)
        except Exception, e:
            logging.warning('Unable to remove file: ' + full_path)
            logging.warning(e)


if __name__ == '__main__':
    # XXX raise KeyError(key)

    parser = argparse.ArgumentParser(prog='do-test')
    parser.add_argument('-r', action='store_true',
                        help='Replace test results')
    parser.add_argument('testname', nargs='*',
                        help='test to run')
    args = parser.parse_args()

    replace_ref = args.r

    # XXX - get directories from environment or command line options

    tester = Tester()

    status = tester.set_environment(
        os.environ['TEST_INPUT_DIR'], 
        os.environ['TEST_OUTPUT_DIR'],
        os.environ['EXECUTABLE_DIR'])
    if not status:
        raise SystemExit('Error in environment')
    tester.open_log('test.log')

    if args.testname:
        test_names = args.testname
    else:
        test_names = [ 'tutorial', 'example', 'include', 'names', 'strings' ]

    logging.info('Tests to run: {}'.format( ' '.join(test_names)))

    pass_names = []
    fail_names = []
    for name in test_names:
        status = tester.set_test(name, replace_ref)

        if status:
            status = tester.do_test()

            if status and not replace_ref:
                status = tester.do_compare()

        if status:
            pass_names.append(name)
            print('{} pass'.format(name))
        else:
            fail_names.append(name)
            print('{} failed'.format(name))

    # summarize results
    if fail_names:
        exit_status = 1
        msg = "Not all tests passed"
    else:
        exit_status = 0
        msg = "All tests passed"
    print(msg)
    logging.info(msg)

    tester.close_log()
    sys.exit(exit_status)
