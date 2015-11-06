#!/bin/bash
###############################################################################
# helper script to change change group and access perms for toolkit
# third party libs
###############################################################################
export LC_TPL_PATH=/usr/gapps/asctoolkit/thirdparty_libs
# set group to toolkitd
chgrp -R toolkitd  ${LC_TPL_PATH}/spack ${LC_TPL_PATH}/*.cmake
# allow group members to read, write, and exec
chmod -R g+rwX ${LC_TPL_PATH}/spack ${LC_TPL_PATH}/*.cmake
# allow everyone else to read, and exec
chmod -R a+rX ${LC_TPL_PATH}/spack ${LC_TPL_PATH}/*.cmake

