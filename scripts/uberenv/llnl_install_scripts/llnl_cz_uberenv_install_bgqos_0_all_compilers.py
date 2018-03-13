#!/usr/local/bin/python

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
 file: llnl_bgq_uberenv_install_bgq_all_compilers.py

 description: 
  uses uberenv to install tpls for the set of compilers we want
  for llnl bgq platforms.

"""

from llnl_lc_uberenv_install_tools import *

def main():
    builds_dir = "/usr/workspace/wsa/axom/thirdparty_libs/builds/"
    specs = ["%clang@4.0.0~cmake~devtools~python~lua"]
    return full_build_and_test_of_tpls(builds_dir,specs)

if __name__ == "__main__":
    sys.exit(main())
