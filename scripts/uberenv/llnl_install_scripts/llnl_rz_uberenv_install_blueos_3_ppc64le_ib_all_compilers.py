#!/usr/local/bin/python
###############################################################################
# 
###############################################################################

"""
 file: llnl_rz_uberenv_blueos_3_ppc64le_ib_all_compilers.py

 description: 
  uses uberenv to install tpls for the set of compilers we want
  for llnl blueos platforms.

"""

from llnl_lc_uberenv_install_tools import *

def main():
    builds_dir = "/usr/workspace/wsrzc/axom/thirdparty_libs/builds/"
    specs = ["%gcc@4.9.3",
             "%clang@coral"]
    return full_build_and_test_of_tpls(builds_dir,specs)

if __name__ == "__main__":
    sys.exit(main())
