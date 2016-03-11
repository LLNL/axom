#!/bin/bash
###############################################################################
# use uberenv to build third party libs on chaos 5 using clang 3.5.0
###############################################################################
source setup_dir.sh
python ../uberenv.py --prefix $LC_TPL_PATH --spec %clang@3.5.0
# set the proper permissions
./lc_install_set_perms.sh

