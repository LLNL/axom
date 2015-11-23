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
import filecmp
import logging
import os
import subprocess
import sys

#subprocess.call("ls -l", shell=True)

#proc = subprocess.Popen(['tail', '-500', 'mylogfile.log'], stdout=subprocess.PIPE)
#for line in proc.stdout.readlines():
#    print line.rstrip()


test_source_dir = ''
test_binary_dir = ''

executable_output_path = ''
code_path = ''

def do_test(name, replace_ref):
    """ Run test, return True/False for pass/fail.
    Files must compare, with no extra or missing files.
    """
    status = True  # assume it passes
    logging.info('--------------------------------------------------')
    logging.info('Testing ' + name)

    testyaml = os.path.join(test_source_dir, name + '.yaml')
    logging.info('Input file: ' + testyaml)
    if not os.path.isfile(testyaml):
        logging.error('Input file does not exist')
        return False

    ref_dir = os.path.join(test_source_dir, name)
    logging.info('Reference directory: ' + ref_dir)

    if replace_ref:
        # replacing reference, just create directly in ref directory
        result_dir = ref_dir
    else:
        result_dir = os.path.join(test_binary_dir, 'tests', name)
        logging.info('Result directory: ' + result_dir)
        makedirs(result_dir)
        clear_files(result_dir)

    output_file = os.path.join(result_dir, 'output')

    cmd = [
        code_path,
        '--path', test_source_dir,
        '--logdir', result_dir,
        '--outdir', result_dir,
        testyaml,
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
    fp = open(output_file, 'w')
    fp.write(output)
    fp.close()

    """
    pipes = subprocess.Popen(cmnd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std_out, std_err = pipes.communicate()

    if pipes.returncode != 0:
        // an error happened!
        err_msg = "%s. Code: %s" % (std_err.strip(), pipes.returncode)
        raise Exception(err_msg)

    elif len(std_err):
        // return code is 0 (no error), but we may want to
        // do something with the info on std_err
        // i.e. logger.warning(std_err)

    // do whatever you want with std_out
    // i.e. json.loads(std_out)
"""

    # collect output file names
#    outfiles = []
#    for line in output.splitlines():
#        if line.startswith('Wrote '):
#            outfiles.append( line[6:] )
#    print("XXXX", outfiles)

    if replace_ref:
        return True

    cmp = filecmp.dircmp(ref_dir, result_dir)
    if not os.path.exists(ref_dir):
        logging.info('Reference directory does not exist: ' + ref_dir)
        return False

    match, mismatch, errors = filecmp.cmpfiles(ref_dir, result_dir, cmp.common)
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
        logging.info('Test {} pass'.format(name))
    else:
        logging.info('Test {} fail'.format(name))
    return status


def makedirs(path):
    """ Make sure directory exists.
    """
    try: 
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise
    # os.makedirs(path,exist_ok=True) python3  3.2


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

    test_binary_dir = os.environ['TEST_BINARY_DIR']
    test_source_dir = os.environ['TEST_SOURCE_DIR']
    executable_output_path = os.environ['EXECUTABLE_OUTPUT_PATH']
    # XXX check existence of directories

    if not os.path.isdir(test_binary_dir):
        raise SystemExit('Missing binary directory: ' + test_binary_dir)
    if not os.path.isdir(test_source_dir):
        raise SystemExit('Missing source directory: ' + test_source_dir)
    if not os.path.isdir(executable_output_path):
        raise SystemExit('Missing executable directory: ' + executable_output_path)

    logname = 'test.log'
    logging.basicConfig(filename=os.path.join(test_binary_dir, logname),
                        filemode='w',
                        level=logging.DEBUG,
                        )

    if args.testname:
        test_names = args.testname
    else:
        test_names = [ 'tutorial', 'example', 'names' ]

    logging.info('Tests to run: {}'.format( ' '.join(test_names)))


    code_path = os.path.join(executable_output_path, 'shroud')
    logging.info('Code to test: ' + code_path)

    result_dir = os.path.join(test_binary_dir, 'tests')
    makedirs(result_dir)

    pass_names = []
    fail_names = []
    for name in test_names:
        status = do_test(name, replace_ref)
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

    logging.shutdown()
    sys.exit(exit_status)
