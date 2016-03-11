#!/bin/bash
###############################################################################
# use uberenv to build third party libs on chaos 5 using gcc 4.9.3
###############################################################################
#use umask to preserve some group perms
source lc_umask.sh
#call uberenv to install our packages
python ../uberenv.py --force --prefix /usr/gapps/asctoolkit/thirdparty_libs/nightly/ --spec %gcc@4.9.3
# change perms as well since umask doesn't help us preserve grwx on newly created dirs
./lc_install_set_perms.sh
