#!/bin/bash
###############################################################################
# use uberenv to build third party libs on bgq using clang 3.7.0
###############################################################################
python ../uberenv.py --prefix /usr/gapps/asctoolkit/thirdparty_libs/ --spec %clang@3.7.0
# set the proper permissions
./lc_install_set_perms.sh


