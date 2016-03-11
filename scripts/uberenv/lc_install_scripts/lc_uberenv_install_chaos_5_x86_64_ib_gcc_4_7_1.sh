#!/bin/bash
###############################################################################
# use uberenv to build third party libs on chaos 5 using gcc 4.7.1
###############################################################################
source setup_dir.sh
python ../uberenv.py --prefix $LC_TPL_PATH --spec %gcc@4.7.1
# set the proper permissions
./lc_install_set_perms.sh

