#!/usr/local/bin/python
###############################################################################
# 
###############################################################################

"""
 file: llnl_lc_find_unused_tpl_dirs.py

 description:

  Finds axom tpl installs on lc that are not referenced in host configs.
  Saves json summary to the file "tpl_dirs_summary.json"

  The json file contains a dictonary with the following entries:

      referenced: contains dirs referenced in host-config files
      found: contains dirs found on the current system
      active: contains dirs both found and referenced 
      unused: contains dirs found but not referenced 

  Note: it is possible to have dirs that are referenced but not
  found (or active, unused) b/c host-configs reference dirs
  across the center (ex: cz + rz)

"""

import glob
import os
import json
import shutil
import sys

from os.path import join as pjoin

from llnl_lc_uberenv_install_tools import sexe, timestamp


def unique(lst):
    """
    Returns list with unique entries of input list
    """
    return list(set(lst))


def lc_tpl_root_dirs():
    """
    Returns list of root dirs axom uses for tpl installs on lc
    """
    return ["/usr/workspace/wsa/axom/thirdparty_libs/builds/", 
            "/usr/workspace/wsrzc/axom/thirdparty_libs/builds/"]    

def clone_axom():
    """
    Creates a fresh clone of axom
    """
    cwd = os.getcwd()
    tmp_dir = os.path.abspath("_tmp_clone_%s" % (timestamp()))
    if os.path.isdir(tmp_dir):
        shutil.rmdir(tmp_dir)
    os.mkdir(tmp_dir)
    os.chdir(tmp_dir)
    print "[cloning axom into %s]" % pjoin(tmp_dir,"axom")
    res = sexe("git clone ssh://git@cz-bitbucket.llnl.gov:7999/atk/axom.git",echo=True)
    if res != 0:
        print "[ERROR: clone of axom repo failed]"
        sys.exit(res)
    os.chdir(cwd)
    return tmp_dir


def list_git_branches():
    """
    Returns a list of all remote git branches
    """
    branches = []
    rcode, branches_out = sexe("git branch -a",
                               ret_output=True)
    for l in branches_out.split():
        if l.startswith("remotes/origin"):
            branch_name = l[l.find("remotes/origin"):]
            branches.append(branch_name)
    return branches

def find_tpl_dirs_in_host_configs_for_all_branches():
    """
    Returns a list of tpl dirs referenced in host configs 
    in all of axom's git branches. 
    """
    res = []
    for branch in list_git_branches():
        print "\n[checking host configs in branch %s]" % branch
        sexe("git checkout %s" % branch,echo=True)
        for tpl_dir in find_tpl_dirs_in_host_configs():
            res.append(tpl_dir)
    return unique(res)
    

def find_tpl_dirs_in_host_configs():
    """
    Returns a list of tpl dirs referenced in host configs 
    in axom.
    """
    prefixes = lc_tpl_root_dirs()
    res = []
    for hc in glob.glob("host-configs/*.cmake"):
        for l in open(hc).readlines():
            for p in prefixes:
               if p in l:
                   ent = l[l.find(p):]
                   ent = ent[:ent.find("spack")]
                   ent = os.path.abspath(ent)
                   res.append(ent)
    return unique(res)

def list_installed_tpl_dirs():
    """
    Returns a list of all lc axom tpl dirs
    """
    prefixes = lc_tpl_root_dirs()
    host_configs = []
    for prefix in prefixes:
        hc_pat = pjoin(prefix,"*/*.cmake")
        for p in glob.glob(hc_pat):
            ent = os.path.split(p)[0]
            ent = os.path.abspath(ent)
            host_configs.append(ent)
    return unique(host_configs)


def check_installed_tpl_dirs():
    """
    Returns a dictionary that contains info about all found tpl 
    dirs and which of these are actually being used by checked-in
    host configs.

      referenced: contains dirs referenced in host-config files
      found: contains dirs found on the current system
      active: contains dirs both found and referenced 
      unused: contains dirs found but not referenced 

    Note: it is possible to have dirs that are referenced but not
    found (or active, unused) b/c host-configs reference dirs
    across the center (ex: cz + rz)
    """
    res = {"referenced":[],
           "found":[],
           "active": [], 
           "unused": []}

    installed_dirs = list_installed_tpl_dirs()
    res["found"].extend(installed_dirs)  

    cwd = os.getcwd()
    tmp_dir = clone_axom()
    os.chdir(pjoin(tmp_dir,"axom"))
    refed_dirs = find_tpl_dirs_in_host_configs_for_all_branches()
    res["referenced"].extend(refed_dirs)

    os.chdir(cwd)
    print "[cleaning up %s]" % tmp_dir
    shutil.rmtree(tmp_dir)

    for l in res["found"]:
        if not l in res["referenced"]:
           res["unused"].append(l)
        else:
           res["active"].append(l)
    return res


def main():
    summary_file ="tpl_dirs_summary.json"
    r = check_installed_tpl_dirs()
    rjson= json.dumps(r,indent=2)
    open(summary_file,"w").write(rjson)
    print rjson
    print ""
    print "[# of referenced %d ]" % len(r["referenced"])
    print "[# of found %d ]" % len(r["found"])
    print "[# of active %d ]" % len(r["active"])
    print "[# of unused %d ]" % len(r["unused"])


if __name__ == "__main__":
    main()









