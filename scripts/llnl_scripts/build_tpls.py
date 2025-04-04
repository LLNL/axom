#!/bin/sh
"exec" "python3" "-u" "-B" "$0" "$@"

# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

"""
 file: build_tpls.py

 description: 
  uses uberenv to install tpls for the set of compilers we want
  for current machine.

"""

from llnl_lc_build_tools import *

from argparse import ArgumentParser

import os


def parse_args():
    "Parses args from command line"
    parser = ArgumentParser()
    # Directory to do all the building
    parser.add_argument("-d", "--directory",
                        dest="directory",
                        default="",
                        help="Location to build all TPL's, timestamp directory will be created (Defaults to shared location)")
    # Spack spec to use for the build
    parser.add_argument("-s", "--spec",
                        dest="spec",
                        default="",
                        help="Spack spec to build (defaults to all available on SYS_TYPE)")
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        dest="verbose",
                        default=False,
                        help="Output logs to screen as well as to files")
    parser.add_argument("-m", "--mirror",
                        dest="mirror",
                        default="",
                        help="Mirror location to use (defaults to shared location)")

    ###############
    # parse args
    ###############
    opts = parser.parse_args()
    # we want a dict b/c the values could 
    # be passed without using argparse
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
        builds_dir = get_shared_libs_dir()
    builds_dir = os.path.abspath(builds_dir)

    repo_dir = get_repo_dir()

    try:
        original_wd = os.getcwd()
        os.chdir(repo_dir)

        timestamp = get_timestamp()
        res = full_build_and_test_of_tpls(builds_dir, timestamp, opts["spec"], opts["verbose"], opts["mirror"])
    finally:
        os.chdir(original_wd)

    return res


if __name__ == "__main__":
    sys.exit(main())
