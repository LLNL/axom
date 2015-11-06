#!/bin/bash
###############################################################################
# helper script to change change group and access perms for toolkit
# third party libs
###############################################################################
# set group to toolkitd
chgrp -R toolkitd /usr/gapps/asctoolkit/thirdparty_libs/
# allow group members to read, write, and exec
chmod -R g+rwX /usr/gapps/asctoolkit/thirdparty_libs/
# allow everyone else to read, and excc
chmod -R a+rX /usr/gapps/asctoolkit/thirdparty_libs/

