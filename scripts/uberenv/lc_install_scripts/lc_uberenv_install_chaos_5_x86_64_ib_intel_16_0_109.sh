#!/bin/bash
###############################################################################
# use uberenv to build third party libs on chaos 5 using gcc 16.0.0
###############################################################################
python ../uberenv.py --prefix $LC_TPL_PATH --spec %intel@16.0.109
# set the proper permissions
./lc_install_set_perms.sh

