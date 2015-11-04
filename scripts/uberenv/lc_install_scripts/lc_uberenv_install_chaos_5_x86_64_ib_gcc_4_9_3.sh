#!/bin/bash
###############################################################################
# use uberenv to build third party libs on chaos 5 using gcc 4.9.3
###############################################################################
python ../uberenv.py --prefix /usr/gapps/asctoolkit/thirdparty_libs/ --spec %gcc@4.9.3
# set the proper permissions
./lc_install_set_perms.sh
