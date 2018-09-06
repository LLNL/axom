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
 file: archive_job.py

 description: 
  Archives results from a non-source or TPL job. Passing an exit code to this
  script will preserve the correct exit status of the automated job

"""

import os
import sys

from llnl_lc_build_tools import *
from optparse import OptionParser


def parse_args():
    "Parses args from command line"
    parser = OptionParser()
    
    parser.add_option("-e", "--exitcode",
                      type="int",
                      dest="exitcode",
                      default="0",
                      help="Exit code to be returned from this script (Defaults to 0)")

    parser.add_option("-n", "--name",
                      dest="name",
                      default="",
                      help="Archive build results under given name")

    parser.add_option("-f", "--files",
                      dest="files",
                      default="",
                      help="Files to be archived (comma delimited)")

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

    timestamp = get_timestamp()
    archive_dir = pjoin(get_archive_base_dir(), get_system_type())
    archive_dir = pjoin(archive_dir, normalize_job_name(opts["name"]), timestamp)

    cwd = os.getcwd()

    print "[Starting Archiving]"
    print "[  Archive Dir: %s]" % archive_dir
    print "[  Job Dir: %s]" % cwd

    # Create a success/failure file if not already created by the job
    if (opts["exitcode"] != 0) and not os.path.exists("failed.json"):
        log_failure(cwd, opts["name"], timestamp)
    elif (opts["exitcode"] == 0) and not os.path.exists("success.json"):
        log_success(cwd, opts["name"], timestamp)

    # Copy any existing info files
    if not os.path.exists(archive_dir):
        os.makedirs(archive_dir)

    copy_if_exists(pjoin(cwd, "info.json"), archive_dir)
    copy_if_exists(pjoin(cwd, "failed.json"), archive_dir)
    copy_if_exists(pjoin(cwd, "success.json"), archive_dir)

    # Copy explicitly given files
    for file in opts["files"].split(','):
        filepath = file.strip()
        if filepath != "":
            copy_if_exists(file, archive_dir)

    return opts["exitcode"]


if __name__ == "__main__":
    # Always exit with the error code provided by user to
    # preserve for reporting the automated job status
    exitcode = main()
    sys.exit(exitcode)
