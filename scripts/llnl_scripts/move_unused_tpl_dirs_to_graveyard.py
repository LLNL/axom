#!/bin/sh
"exec" "python" "-u" "-B" "$0" "$@"

###############################################################################
# Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-741217
#
# All rights reserved.
#
# This file is part of Axom.
#
# For details about use and distribution, please read axom/LICENSE.
###############################################################################


"""
 file: llnl_lc_move_unused_tpl_dirs_to_graveyard.py

 description:

  Reads JSON result from llnl_lc_find_unused_tpl_dirs.py,
  and moves those TPL dirs marked "unused" to the 
  "_graveyard" for cleanup.

"""

import os
import json

from os.path import join as pjoin

from llnl_lc_build_tools import sexe


def graveyard_path(d):
    # TODO harden +s test on RZ
    return pjoin(d[:d.rfind("builds/")+7],"_graveyard")

def main():
    summary_file ="tpl_dirs_summary.json"
    r= json.load(open(summary_file))
    print ""
    print "[# of referenced %d ]" % len(r["referenced"])
    print "[# of found %d ]" % len(r["found"])
    print "[# of active %d ]" % len(r["active"])
    print "[# of unused %d ]" % len(r["unused"])
    for d in r["unused"]:
        cmd = "mv %s %s" % (d, graveyard_path(d)) 
        sexe(cmd)


if __name__ == "__main__":
    main()









