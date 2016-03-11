#!/bin/bash
###############################################################################
# use uberenv to build third party libs on chaos 5 using gcc 4.9.3
###############################################################################
# create new library directory
source setup_dir.sh
python ../uberenv.py --prefix $LC_TPL_PATH --spec %gcc@4.9.3
#set the proper permissions
./lc_install_set_perms.sh
