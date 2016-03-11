#!/bin/bash
###############################################################################
# helper script to change change group and access perms for toolkit
# third party libs
###############################################################################
# set group to toolkitd
echo "Changing group to toolkit for $LC_TPL_PATH/spack, $LC_TPL_PATH/*.cmake"
chgrp -f -R toolkitd  ${LC_TPL_PATH}/spack ${LC_TPL_PATH}/*.cmake
# allow group members to read, write, and exec
echo "Allowing group read, write, execute on $LC_TPL_PATH/spack, $LC_TPL_PATH/*.cmake"
chmod -f -R g+rwX ${LC_TPL_PATH}/spack ${LC_TPL_PATH}/*.cmake
# allow everyone else to read, and exec
echo "Allowing everyone read and execute on $LC_TPL_PATH/spack, $LC_TPL_PATH/*.cmake"
chmod -f -R a+rX ${LC_TPL_PATH}/spack ${LC_TPL_PATH}/*.cmake

