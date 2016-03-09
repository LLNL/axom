#!/bin/bash
###############################################################################
# use uberenv to build third party libs on bgq using clang 3.7.0
###############################################################################
#use umask to preserve group perms
source lc_umask.sh
#call uberenv to install our packages
python ../uberenv.py --prefix /usr/gapps/asctoolkit/thirdparty_libs/nightly --spec %clang@3.7.0


