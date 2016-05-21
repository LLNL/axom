#!/usr/l:qocal/bin/python
###############################################################################
# 
###############################################################################

"""
 file: llnl_lc_uberenv_test_host_configs.py

 description: 
  builds and tests toolkit against a set of host config files created by
  spack / uberenv. 

  To use, set the env var UBERENV_PREFIX to the path of an uberenv install.

"""

import os
import sys

from llnl_lc_uberenv_install_tools import *


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





