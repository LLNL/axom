#!/bin/bash
###############################################################################
# helper script to change change group and access perms for toolkit
# third party libs
###############################################################################
export LC_TPL_PATH=/usr/gapps/asctoolkit/thirdparty_libs
echo "[changing permissions for $LC_TPL_PATH]"
# set group to toolkitd
chgrp -f -R toolkitd  ${LC_TPL_PATH}/spack ${LC_TPL_PATH}/*.cmake
# allow group members to read, write, and exec
chmod -f -R g+rwX ${LC_TPL_PATH}/spack ${LC_TPL_PATH}/*.cmake
# allow everyone else to read, and exec
chmod -f -R a+rX ${LC_TPL_PATH}/spack ${LC_TPL_PATH}/*.cmake
echo "[finished changing permissions for $LC_TPL_PATH]"

