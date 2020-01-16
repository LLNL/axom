#!/bin/sh
"exec" "python" "-u" "-B" "$0" "$@"

# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

"""
 file: fix_permissions.py

 description: 
  Fixes the file permissions for archiving and tpl locations

"""

import sys
from llnl_lc_build_tools import *

def main():
    for directory in [get_archive_base_dir(), get_shared_tpl_base_dir()]:
        print "Fixing group to axomdev for {0}".format(directory)
        sexe("chgrp -R axomdev {0}".format(directory))
        print "Fixing file permissions for {0}".format(directory)
        sexe("chmod -R g+rwX {0}".format(directory))

    return True


if __name__ == "__main__":
    if main():
        sys.exit(0)
    sys.exit(1)
