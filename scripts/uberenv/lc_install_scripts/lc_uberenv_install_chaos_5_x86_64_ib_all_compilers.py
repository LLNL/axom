###############################################################################
# 
###############################################################################

"""
 file: lc_uberenv_install_chaos5_x86_64_ib_all_compilers.py

 description: 
  uses uberenv to install tpls for the set of compilers we want
  for lc chaos 5 platforms.

"""

import subprocess
import datetime

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

def set_toolkit_fs_perms(prefix):
    # change group to  toolktid
    print "[change grou for toolkitd members to rwX]"
    sexe("chgrp -f -R toolkitd  %s" %(prefix),echo=True)
    print "[change perms for toolkitd members to rwX]"
    # change group perms to rwX
    sexe("chmod -f -R g+rwX  %s" %(prefix),echo=True)
    # change perms for all to rX
    print "[change perms for all users to rX]"
    sexe("chmod -f -R g+rwX  %s" %(prefix),echo=True)

def install_spec(prefix,spec):
    cmd = "python ../uberenv.py --prefix %s --spec %s " % (prefix,spec)
    return sexe(cmd,echo=True)


def main():
    prefix = "/usr/gapps/asctoolkit/thirdparty_libs/builds/" + timestamp()
    specs = ["%clang@3.5.0",
             "%gcc@4.7.1",
             "%gcc@4.9.3",
             "%intel@15.0.187",
             "%intel@16.0.109"]
    for spec in specs:
        install_spec(prefix,spec)
    set_toolkit_fs_perms(prefix)


if __name__ == "__main__":
    main()


