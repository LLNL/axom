#!/usr/l:qocal/bin/python
###############################################################################
# 
###############################################################################

"""
 file: llnl_lc_uberenv_test_host_configs.py

 description: 
  builds and tests toolkit against a set of host config files created by
  spack / uberenv. 

  To use
  UBERENV_PREFIX

"""

import os
import sys

from llnl_lc_uberenv_install_tools import *

############################################################
# helpers for testing a set of host configs
############################################################

def build_and_test_host_config(test_root,host_config):
    host_config_root = os.path.splitext(os.path.basename(host_config))[0]
    # setup build and install dirs
    build_dir   = pjoin(test_root,"build-%s"   % host_config_root)
    install_dir = pjoin(test_root,"install-%s" % host_config_root)
    # configure
    sexe("python ../../config-build.py  -bp %s -ip %s -hc %s" % (build_dir,install_dir,host_config),echo=True)
    # build
    sexe("cd %s && make -j 8" % build_dir,echo=True)
    sexe("cd %s && make test" % build_dir,echo=True)
    sexe("cd %s && make install" % build_dir,echo=True)
    # check for results in make install?


def build_and_test_host_configs(prefix):
    host_configs = glob.glob(pjoin(prefix,"*.cmake"))
    if len(host_configs) > 0:
        test_root =  pjoin(prefix,"_asctk_build_and_test_%s" % timestamp())
        os.mkdir(test_root)
        write_build_info(pjoin(test_root,"info.json")) 
        for host_config in host_configs:
            build_and_test_host_config(test_root,host_config)
        set_toolkit_group_and_perms(test_root)
    else:
        print "[error no host configs found at %s]" % prefix


def usage():
    print "[usage: env UBERENV_PREFIX={path} python build_and_test_host_configs.py]"

def main():
    if not "UBERENV_PREFIX" in os.environ.keys():
        usage()
        sys.exit(-1)
    prefix = os.environ["UBERENV_PREFIX"]
    prefix = os.path.abspath(prefix)
    if not os.path.isdir(prefix):
        print "[error %s is not a valid directory]" % prefix
        sys.exit(-1)
    build_and_test_host_configs(prefix)


if __name__ == "__main__":
    main()





