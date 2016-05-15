#!/usr/local/bin/python
###############################################################################
# 
###############################################################################

"""
 file: llnl_lc_set_toolkit_permissions.py

 description: 
  sets group and user perms for toolkit on llnl lc systems.

"""

import subprocess

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

def set_toolkit_perms(fdir):
    print "[changing group and access perms of: %s]" % fdir
    # change group to toolktid
    print "[changing group to toolkitd]"
    sexe("chgrp -f -R toolkitd  %s" % (fdir),echo=True)
    # change group perms to rwX
    print "[changing perms for toolkitd members to rwX]"
    sexe("chmod -f -R g+rwX  %s" % (fdir),echo=True)
    # change perms for all to rX
    print "[changing perms for all users to rX]"
    sexe("chmod -f -R a+rX  %s" % (fdir),echo=True)


def main():
    if len(sys.argv) <1:
        print "usage: llnl_lc_set_toolkit_perms.py [directory]"
        sys.exit(-1)
    fdir = sys.argv[1]
    set_toolkit_perms(fdir)

if __name__ == "__main__":
    main()


