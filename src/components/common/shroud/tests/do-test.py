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

from __future__ import print_function

import os
import argparse
import filecmp
import subprocess

#subprocess.call("ls -l", shell=True)

#proc = subprocess.Popen(['tail', '-500', 'mylogfile.log'], stdout=subprocess.PIPE)
#for line in proc.stdout.readlines():
#    print line.rstrip()


test_source_dir = ''
test_binary_dir = ''

executable_output_path = ''
code_path = ''

def do_test(name):
    # XXX make sure input file exists
    ref_dir = os.path.join(test_source_dir, name)
    # make sure output directory exists
    result_dir = os.path.join(test_binary_dir, 'tests', name)
    makedirs(result_dir)
    # XXX remove old generated files

    cmd = [
        code_path,
        '--logdir', result_dir,
        '--outdir', result_dir,
        os.path.join(test_source_dir, name + '.yaml'),
        ]
    print("LOG", ' '.join(cmd))


    try:
        output = subprocess.check_output(
            cmd,
            stderr=subprocess.STDOUT,
            universal_newlines=True)
    except subprocess.CalledProcessError as exc:
        print("Status : FAIL", exc.returncode, exc.output)
        raise SystemExit("shroud failed")
#    except subprocess.CalledProcessError:
#        print(output)


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

    print("Comparing", ref_dir, result_dir)
    cmp = filecmp.dircmp(ref_dir, result_dir)
    print("COMMON", cmp.common)
    match, mismatch, errors = filecmp.cmpfiles(ref_dir, result_dir, cmp.common)
    print("MATCH", match)
    print("MISMATCH", mismatch)
    print("ERRORS", errors)
    

        

    print("REF_ONLY", cmp.left_only)
    print("RESULT_ONLY", cmp.right_only)

    raise SystemExit




def call_command(command):
    process = subprocess.Popen(command.split(' '),
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    return process.communicate()

def call_commandlst(command):
    print("call_commandlst", command)
    process = subprocess.Popen(command,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    return process.communicate()


def makedirs(path):
    try: 
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise
    # os.makedirs(path,exist_ok=True) python3  3.2


if __name__ == '__main__':
    # XXX raise KeyError(key)
    test_binary_dir = os.environ['TEST_BINARY_DIR']
    test_source_dir = os.environ['TEST_SOURCE_DIR']
    executable_output_path = os.environ['EXECUTABLE_OUTPUT_PATH']
    # XXX check existence of directories

    code_path = os.path.join(executable_output_path, 'shroud')

    result_dir = os.path.join(test_binary_dir, 'tests')
    makedirs(result_dir)
    
    do_test('example')
