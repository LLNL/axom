#!/bin/sh
"exec" "python" "-u" "-B" "$0" "$@"

# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)


"""
 file: move_unused_tpl_dirs_to_graveyard.py

 description:

  Reads JSON result from find_unused_tpl_dirs.py,
  and moves those TPL dirs marked "unused" to the 
  "_graveyard" for cleanup.

"""

import os
import json

from os.path import join as pjoin

from llnl_lc_build_tools import *

gy_base = pjoin(get_shared_libs_dir(), "_graveyard")

def graveyard_path(d):
    path = d.replace(get_shared_libs_dir(), gy_base)
    return path

def main():
    # Create directory structure if it's not there
    for sys_type in get_supported_sys_types():
        directory = pjoin(gy_base, sys_type)
        if not os.path.exists(directory):
            os.makedirs(directory)

    summary_file ="tpl_dirs_summary.json"
    if not os.path.exists(summary_file):
        print("Error: Summary file from find_unused_tpl_dirs.py did not exist: {0}", summary_file)
        return False

    r = json.load(open(summary_file))

    print ""
    print "[# of referenced %d ]" % len(r["referenced"])
    print "[# of found %d ]" % len(r["found"])
    print "[# of active %d ]" % len(r["active"])
    print "[# of unused %d ]" % len(r["unused"])

    success = True
    for d in r["unused"]:
        cmd = "mv %s %s" % (d, graveyard_path(d)) 
        res = sexe(cmd)
        if res != 0:
            success = False
    return success


if __name__ == "__main__":
    if main():
        sys.exit(0)
    sys.exit(1)
