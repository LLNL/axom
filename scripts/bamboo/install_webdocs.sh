#!/bin/bash

#######################################
#  Variables
#######################################

BUILD_TYPE="Debug"
BUILD_PATH="axom_build"
INSTALL_PATH="axom_install"

TOOLKIT_WEB_ROOT="/usr/global/web-pages/lc/www/axom"
DOCS_DIR_OLD="${TOOLKIT_WEB_ROOT}/docs_old"
DOCS_DIR="${TOOLKIT_WEB_ROOT}/docs"

if [[ $HOSTNAME == rz* ]]; then
    echo "Error: this job is to be run on the CZ"
    exit 1
else
    HOST_CONFIGURATION="host-configs/quartz-toss_3_x86_64_ib-gcc@7.3.0.cmake"
fi


#######################################
#  Configure
#######################################

echo "------------------------------Configure Build------------------------------"
./config-build.py -ecc -hc $HOST_CONFIGURATION -bt $BUILD_TYPE -bp $BUILD_PATH -ip $INSTALL_PATH
if [ $? -ne 0 ]; then
    echo "Error: config-build.py failed"
    exit 1
fi
echo "-----------------------------------------------------------------------"


#######################################
#  Build documentation
#######################################

echo "------------------------------Making docs------------------------------"
cd $BUILD_PATH
make VERBOSE=1 -j16 docs
if [ $? -ne 0 ]; then
    echo "Error: 'make docs' failed"
    exit 1
fi
echo "-----------------------------------------------------------------------"


#######################################
#  make install
#######################################

echo "--------------------------make install--------------------------------"
make VERBOSE=1 -j16 install
if [ $? -ne 0 ]; then
    echo "Error: 'make install' failed"
    exit 1
fi
echo "-----------------------------------------------------------------------"


#######################################
#  make install
#######################################

cd ..
echo "---------------------------Installing docs to public-------------------"
if [ -d  ${DOCS_DIR_OLD} ]; then
   rm -rf ${DOCS_DIR_OLD}
   if [ $? -ne 0 ]; then
       echo "Warning: 'rm' of docs_old failed"
   fi
fi

if [ -d  ${DOCS_DIR} ]; then
    mv ${DOCS_DIR} ${DOCS_DIR_OLD}
    if [ $? -ne 0 ]; then
        echo "Warning: 'mv' docs to docs_old failed"
    fi
fi

cp -R ./${INSTALL_PATH}/docs ${TOOLKIT_WEB_ROOT}/
if [ $? -ne 0 ]; then
    echo "Error: 'cp' failed"
    exit 1
fi

chgrp -R axom ${DOCS_DIR}
if [ $? -ne 0 ]; then
    echo "Error: 'chgrp' failed"
    exit 1
fi

chmod -R g+r+w+X ${DOCS_DIR}
if [ $? -ne 0 ]; then
    echo "Error: 'chmod' failed"
    exit 1
fi

chmod -R o+r+X ${DOCS_DIR}
if [ $? -ne 0 ]; then
    echo "Error: 'chmod' failed"
    exit 1
fi
echo "-----------------------------------------------------------------------"
