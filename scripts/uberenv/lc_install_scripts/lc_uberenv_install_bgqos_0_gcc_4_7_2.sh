#!/bin/bash
###############################################################################
# use uberenv to build third party libs on bgq using clang 4.7.2
###############################################################################
python ../uberenv.py --prefix /usr/gapps/asctoolkit/thirdparty_libs/ --spec %gcc@4.7.2
# set the proper permissions
./lc_install_set_perms.sh

