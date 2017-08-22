#!/usr/local/bin/python
###############################################################################
# 
###############################################################################

"""
 file: llnl_rz_uberenv_install_bgqos_0_all_compilers.py

 description: 
  uses uberenv to install tpls for the set of compilers we want
  for llnl rz bgq  platforms.

"""

from llnl_lc_uberenv_install_tools import *

def main():
    builds_dir = "/usr/workspace/wsrzc/axom/thirdparty_libs/builds/"
    specs = ["%gcc@4.7.2~cmake~devtools~python~lua",
             "%clang@3.7.0~cmake~devtools~python~lua"]
    return full_build_and_test_of_tpls(builds_dir,specs)

if __name__ == "__main__":
    sys.exit(main())


