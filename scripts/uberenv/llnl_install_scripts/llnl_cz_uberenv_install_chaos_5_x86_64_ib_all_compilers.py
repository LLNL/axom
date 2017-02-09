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
    builds_dir = "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/"
    mirror_dir = pjoin(builds_dir,"mirror")
    # unique install location
    prefix =  pjoin(builds_dir,timestamp())
    # create a mirror
    uberenv_create_mirror(prefix,mirror_dir)
    # write info about this build
    write_build_info(pjoin(prefix,"info.json"))
    # spack specs for the cz chaos systems
    specs = ["%clang@3.5.0",
             "%gcc@4.7.1",
             "%gcc@4.9.3",
             "%intel@15.0.187",
             "%intel@16.0.109"]
    # use uberenv to install for all specs
    for spec in specs:
        res = uberenv_install_tpls(prefix,spec,mirror_dir)
        if res != 0:
            print "[ERROR: Failed build of tpls for spec %s]" % spec
            return res
        else:
            print "[SUCCESS: Finished build tpls for spec %s]" % spec
    # patch manual edits into host config files
    patch_host_configs(prefix)
    # build the toolkit against the new tpls
    build_and_test_host_configs(prefix)
    # set proper perms for installed tpls
    set_toolkit_group_and_perms(prefix)
    # set proper perms for the mirror files
    set_toolkit_group_and_perms(mirror_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
