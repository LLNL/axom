#!/bin/bash
###############################################################################
# use uberenv to build third party libs on chaos 5 using intel 15
###############################################################################
#use umask to preserve group perms
source lc_umask.sh
#call uberenv to install our packages
python ../uberenv.py --prefix /usr/gapps/asctoolkit/thirdparty_libs/nightly --spec %intel@15.0.0
