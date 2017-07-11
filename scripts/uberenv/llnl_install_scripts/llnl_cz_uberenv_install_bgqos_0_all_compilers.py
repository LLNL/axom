#!/usr/local/bin/python
###############################################################################
# 
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
    specs = ["%gcc@4.7.2~cmake~devtools~python~lua",
             "%clang~cmake~devtools~python~lua"]

if __name__ == "__main__":
    sys.exit(main())
