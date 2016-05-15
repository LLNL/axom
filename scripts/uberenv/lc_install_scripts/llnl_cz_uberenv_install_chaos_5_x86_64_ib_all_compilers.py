#!/usr/local/bin/python
###############################################################################
# 
###############################################################################

"""
 file: llnl_cz_uberenv_install_chaos5_x86_64_ib_all_compilers.py

 description: 
  uses uberenv to install tpls for the set of compilers we want
  for llnl cz chaos 5 platforms.

"""

from llnl_lc_uberenv_install_tools import *

def main():
    # install loc for cz
    prefix = "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/" + timestamp()
    # spack specs for the cz chaos systems
    specs = ["%clang@3.5.0",
             "%gcc@4.7.1",
             "%gcc@4.9.3",
             "%intel@15.0.187",
             "%intel@16.0.109"]
    # use uberenv to install for all specs
    for spec in specs:
        uberenv_install_tpls(prefix,spec)
    # set proper perms for installed tpls
    set_toolkit_group_and_perms(prefix)
    # look for host config files as a sanity check
    check_host_configs(prefix)


if __name__ == "__main__":
    main()
