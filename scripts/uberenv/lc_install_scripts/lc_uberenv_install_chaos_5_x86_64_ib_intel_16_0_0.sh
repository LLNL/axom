#!/bin/bash
###############################################################################
# use uberenv to build third party libs on chaos 5 using intel 16
###############################################################################
python ../uberenv.py --prefix /usr/gapps/asctoolkit/thirdparty_libs/nightly --spec %intel@16.0.0
# set the proper permissions
./lc_install_set_perms.sh

