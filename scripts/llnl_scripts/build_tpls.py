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
    if on_rz():
        builds_dir = "/usr/workspace/wsrzc/axom/thirdparty_libs/builds/"
    else:
        builds_dir = "/usr/workspace/wsa/axom/thirdparty_libs/builds/"

    specs = get_specs_for_current_machine()
    
    return full_build_and_test_of_tpls(builds_dir, specs)

if __name__ == "__main__":
    sys.exit(main())
