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
 file: llnl_rz_uberenv_install_toss_3_x86_64_ib_all_compilers.py

 description: 
  uses uberenv to install tpls for the set of compilers we want
  for llnl rz toss 3 platforms.

"""

from llnl_lc_uberenv_install_tools import *

def main():
    script_dir = os.path.dirname(os.path.realpath(__file__))
    repo_dir = os.path.abspath(os.path.join(script_dir, ".."))

    host_configs = get_host_configs_for_current_machine()
    failed = []
    succeeded = []
    for host_config in host_configs:
        if build_and_test_host_config(repo_dir, host_config) != 0:
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
