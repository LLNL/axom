"""
 file: uberenv.py

 description: uses spack to install the external third party libs used by the asctoolkit.

"""

import os
import sys
import subprocess
import shutil
import socket
import platform

from os import environ as env

from optparse import OptionParser
from os.path import join as pjoin


def sexe(cmd):
    "Basic shell call helper."
    print "[sexe:%s ]" % cmd
    subprocess.call(cmd,shell=True)

def parse_args():
    "Parses args from command line"
    parser = OptionParser()
    #
    parser.add_option("--prefix",
                      dest="prefix",
                      default="uberenv_libs",
                      help="destination dir")
    # parse args
    opts, extras = parser.parse_args()
    # we want a dict b/c the values could 
    # be passed without using optparse
    opts = vars(opts)
    return opts, extras

def main():
    """
    clones and runs spack to setup our third_party libs and
    creates a host-config.cmake file that can be used by the asctoolkit.
    """
    if "darwin" in platform.system().lower():
        dep_tgt = platform.mac_ver()[0]
        dep_tgt = dep_tgt[:dep_tgt.rfind(".")]
        print "[setting MACOSX_DEPLOYMENT_TARGET to %s]" % dep_tgt
        env["MACOSX_DEPLOYMENT_TARGET"] = dep_tgt
    # parse args from command line
    opts, extras = parse_args()
    # get the current working path, and the glob used to identify the 
    # package files we want to hot-copy to spack
    cwd = os.path.abspath(os.getcwd())
    pkgs = pjoin(cwd,"packages","*")
    # setup destination paths
    dest_dir = os.path.abspath(opts["prefix"])
    dest_spack = pjoin(dest_dir,"spack")
    dest_spack_pkgs = pjoin(dest_spack,"var","spack","packages")
    # print a warning if the dest path already exists
    if os.path.isdir(dest_spack):
        print "[warning: destination '%s' already exists]"  % dest_spack
    else:
        os.mkdir(dest_dir)
    # clone spack into the dest path
    os.chdir(dest_dir)
    sexe("git clone https://github.com/scalability-llnl/spack.git")
    # hot-copy our packages into spack
    sexe("cp -Rf %s %s" % (pkgs,dest_spack_pkgs))
    # use the uberenv package to trigger the right builds and build an host-config.cmake file
    sexe("spack/bin/spack install uberenv")
    host_cfg = "%s.cmake" % socket.gethostname()
    if os.path.isfile(host_cfg):
        print "[result host-config file: %s]" % os.path.abspath(host_cfg)


if __name__ == "__main__":
    main()


