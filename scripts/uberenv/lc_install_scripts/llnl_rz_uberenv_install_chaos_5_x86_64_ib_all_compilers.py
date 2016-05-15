#!/usr/local/bin/python
###############################################################################
# 
###############################################################################

"""
 file: llnl_cz_uberenv_install_chaos5_x86_64_ib_all_compilers.py

 description: 
  uses uberenv to install tpls for the set of compilers we want
  for llnl rz chaos 5 platforms.

"""

import subprocess
import datetime

from llnl_lc_set_toolkit_perms import set_toolkit_perms

def sexe(cmd,ret_output=False,echo = False):
    """ Helper for executing shell commands. """
    if echo:
        print "[exe: %s]" % cmd
    if ret_output:
        p = subprocess.Popen(cmd,
                             shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
        res =p.communicate()[0]
        return p.returncode,res
    else:
        return subprocess.call(cmd,shell=True)

def timestamp(t=None,sep="_"):
    """ Creates a timestamp that can easily be included in a filename. """
    if t is None:
        t = datetime.datetime.now()
    sargs = (t.year,t.month,t.day,t.hour,t.minute,t.second)
    sbase = "".join(["%04d",sep,"%02d",sep,"%02d",sep,"%02d",sep,"%02d",sep,"%02d"])
    return  sbase % sargs


def install_tpls(prefix,spec):
    cmd = "python ../uberenv.py --prefix %s --spec %s " % (prefix,spec)
    return sexe(cmd,echo=True)


def main():
    prefix = "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/" + timestamp()
    specs = ["%clang@3.5.0",
             "%gcc@4.7.1",
             "%gcc@4.9.3",
             "%intel@15.0.187",
             "%intel@16.0.109"]
    for spec in specs:
        install_tpls(prefix,spec)
    set_toolkit_perms(prefix)


if __name__ == "__main__":
    main()


