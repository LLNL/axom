#!/bin/sh
"exec" "python" "-u" "-B" "$0" "$@"

###############################################################################
# Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-741217
#
# All rights reserved.
#
# This file is part of Axom.
#
# For details about use and distribution, please read axom/LICENSE.
###############################################################################

"""
 file: build_tpls.py

 description: 
  uses uberenv to install tpls for the set of compilers we want
  for current machine.

"""

from llnl_lc_build_tools import *

from optparse import OptionParser

import os


def parse_args():
    "Parses args from command line"
    parser = OptionParser()
    # Directory to do all the building
    parser.add_option("-d", "--directory",
                      dest="directory",
                      default="",
                      help="Location to build all TPL's, timestamp directory will be created (Defaults to shared location)")
    # Whether to archive results
    parser.add_option("-a", "--archive",
                      dest="archive",
                      default="",
                      help="Archive build results under given name (Defaults to off)")

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

    # Determine location to do all the building
    if opts["directory"] != "":
        builds_dir = opts["directory"]
        if not os.path.exists(builds_dir):
            os.makedirs(builds_dir)
    else:
        builds_dir = get_shared_tpl_builds_dir()
    builds_dir = os.path.abspath(builds_dir)

    # Make sure builds_dir is not located inside of the repository (CMake hates this)
    repo_dir = get_repo_dir()
    if builds_dir.startswith(repo_dir):
        print "Error: TPL build directory cannot be inside of the AXOM repository."
        print "   CMake will not allow you to use any TPL's include directory inside the repository."
        return False

    
    if opts["archive"] != "":
        job_name = opts["archive"]
    else:
        job_name = get_username() + "/" + os.path.basename(__file__)

    try:
        original_wd = os.getcwd()
        os.chdir(repo_dir)

        timestamp = get_timestamp()
        res = full_build_and_test_of_tpls(builds_dir, job_name, timestamp)

        if opts["archive"] != "":
            # Get information for archiving
            archive_base_dir = get_archive_base_dir()
            if opts["archive"] != "":
                job_name = opts["archive"]
            else:
                job_name = get_username() + "/" + os.path.basename(__file__)

            archive_tpl_logs(builds_dir, job_name, timestamp)
    finally:
        os.chdir(original_wd)

    return res


if __name__ == "__main__":
    sys.exit(main())
