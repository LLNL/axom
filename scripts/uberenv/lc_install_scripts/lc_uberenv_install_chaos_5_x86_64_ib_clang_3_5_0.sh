#!/bin/bash
###############################################################################
# use uberenv to build third party libs on chaos 5 using clang 3.5.0
###############################################################################
python ../uberenv.py --prefix /usr/gapps/asctoolkit/thirdparty_libs/nightly --spec %clang@3.5.0
# set the proper permissions
./lc_install_set_perms.sh

