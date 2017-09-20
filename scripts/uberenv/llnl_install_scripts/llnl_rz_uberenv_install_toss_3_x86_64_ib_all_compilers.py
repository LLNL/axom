#!/usr/local/bin/python
###############################################################################
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
    builds_dir = "/usr/workspace/wsrzc/axom/thirdparty_libs/builds/"
    specs = ["%clang@3.9.0",
             "%gcc@4.9.3",
             "%intel@16.0.4",
             "%intel@17.0.0"]
    return full_build_and_test_of_tpls(builds_dir,specs)

if __name__ == "__main__":
    sys.exit(main())
