#!/bin/sh
"exec" "python" "-u" "-B" "$0" "$@"

# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

"""
 file: build_src.py

 description: 
  Builds all Axom with the host-configs for the current machine.

"""

from llnl_lc_build_tools import *

from optparse import OptionParser


def parse_args():
    "Parses args from command line"
    parser = OptionParser()
    # Location of source directory to build
    parser.add_option("-d", "--directory",
                      dest="directory",
                      default="",
                      help="Directory of source to be built (Defaults to current)")
    # Whether to archive results
    parser.add_option("-a", "--archive",
                      dest="archive",
                      default="",
                      help="Archive build results under given name (Defaults to off)")

    # Specify build type
    parser.add_option("-b", "--build_type",
                      dest="build_type",
                      default="Debug",
                      help="Specify build type (Defaults to Debug)")

    ###############
    # parse args
    ###############
    opts, extras = parser.parse_args()
    # we want a dict b/c the values could 
    # be passed without using optparse
    opts = vars(opts)
    return opts


def main():
    opts = parse_args()

    # Determine source directory to be built
    if os.environ.get("UBERENV_PREFIX") != None:
        repo_dir = os.environ["UBERENV_PREFIX"]
        if not os.path.isdir(repo_dir):
            print "[ERROR: Given environment variable 'UBERENV_PREFIX' is not a valid directory]"
            print "[    'UBERENV_PREFIX' = %s]" % repo_dir
            return 1
    if opts["directory"] != "":
        repo_dir = opts["directory"]
        if not os.path.isdir(repo_dir):
            print "[ERROR: Given command line variable '--directory' is not a valid directory]"
            print "[    '--directory' = %s]" % repo_dir
            return 1
    else:
        repo_dir = get_repo_dir()

    if opts["archive"] != "":
        job_name = opts["archive"]
    else:
        job_name = get_username() + "/" + os.path.basename(__file__)

    if opts["build_type"] in ["Release", "Debug", "RelWithDebInfo", "MinSizeRel"]:
        build_type = opts["build_type"]
    else:
        build_type= "Debug"

    try:
        original_wd = os.getcwd()
        os.chdir(repo_dir)
        timestamp = get_timestamp()
        res = build_and_test_host_configs(repo_dir, job_name, build_type, timestamp, False)

        # Archive logs
        if opts["archive"] != "":
            archive_src_logs(repo_dir, job_name, timestamp)
    finally:
        os.chdir(original_wd)

    return res

if __name__ == "__main__":
    sys.exit(main())
