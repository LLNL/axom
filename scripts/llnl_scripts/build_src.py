#!/bin/sh
"exec" "python" "-u" "-B" "$0" "$@"

###############################################################################
# Copyright (c) 2017, Lawrence Livermore National Security, LLC.
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
# 
###############################################################################

"""
 file: build_src.py

 description: 
  Builds all Axom with the host-configs for the current machine.

"""

from llnl_lc_uberenv_install_tools import *

from optparse import OptionParser


def parse_args():
    "Parses args from command line"
    parser = OptionParser()
    # where to install
    parser.add_option("--src",
                      dest="src",
                      default="",
                      help="Source to be built (Defaults to current)")

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

    if os.environ["UBERENV_PREFIX"] != "":
        src_dir = os.environ["UBERENV_PREFIX"]
        if not os.path.isdir(src_dir):
            print "[ERROR: Given environment variable 'UBERENV_PREFIX' is not a valid directory]"
            print "[    'UBERENV_PREFIX' = %s" % src_dir
            return 1
    if opts["src"] != "":
        src_dir = opts["src"]
        if not os.path.isdir(src_dir):
            print "[ERROR: Given command line variable '--src' is not a valid directory]"
            print "[    '--src' = %s" % src_dir
            return 1
    else:
        script_dir = os.path.dirname(os.path.realpath(__file__))
        src_dir = os.path.abspath(os.path.join(script_dir, "../.."))

    host_configs = get_host_configs_for_current_machine(src_dir)
    failed = []
    succeeded = []
    for host_config in host_configs:
        if build_and_test_host_config(src_dir, host_config) != 0:
            failed.append(host_config)
        else:
            succeeded.append(host_config)

    print "\n\n\n"

    if len(succeeded) > 0:
        print "Succeeded:"
        for host_config in succeeded:
            print "    " + host_config

    if len(failed) > 0:
        print "Failed:"
        for host_config in failed:
            print "    " + host_config
        print "\n"
        return 1

    print "\n"
    return 0

if __name__ == "__main__":
    sys.exit(main())
