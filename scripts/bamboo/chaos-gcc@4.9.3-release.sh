#!/bin/bash


COMPILER="gcc@4.9.3"

if [[ $HOSTNAME == rz* ]]; then
    HC="host-configs/rzmerl-chaos_5_x86_64_ib-${COMPILER}.cmake"
else
    HC="host-configs/surface-chaos_5_x86_64_ib-${COMPILER}.cmake"
fi

BT="Release"
#Note (KW): for some reason, ${COMPILER} was not present in the BP or IP here, 
#           but was present in all other scripts. Should we add it?
BP="build-chaos-${BT,,}"
IP="install-chaos-${BT,,}"
COMP_OPT=""
BUILD_OPT=""


echo "Configuring..."
echo "-----------------------------------------------------------------------"
./scripts/config-build.py -ecc -hc $HC -bt $BT -bp $BP -ip $IP $COMP_OPT $BUILD_OPT    
if [ $? -ne 0 ]; then
    echo "Error: config-build.py failed"
    exit 1
fi
echo "-----------------------------------------------------------------------"

cd $BP
    echo "Building..."
    echo "-----------------------------------------------------------------------"
    make VERBOSE=1 -j16
    if [ $? -ne 0 ]; then
        echo "Error: 'make' failed"
        exit 1
    fi
    echo "-----------------------------------------------------------------------"

    echo "Running tests..."
    echo "-----------------------------------------------------------------------"
    make test ARGS="-T Test"
    if [ $? -ne 0 ]; then
        echo "Error: 'make test' failed"
        exit 1
    fi
    echo "-----------------------------------------------------------------------"

    echo "Making docs..."
    echo "-----------------------------------------------------------------------"
    make VERBOSE=1 docs
    if [ $? -ne 0 ]; then
        echo "Error: 'make docs' failed"
        exit 1
    fi
    echo "-----------------------------------------------------------------------"

    echo "Installing files..."
    echo "-----------------------------------------------------------------------"
    make VERBOSE=1 install
    if [ $? -ne 0 ]; then
        echo "Error: 'make install' failed"
        exit 1
    fi
    echo "-----------------------------------------------------------------------"
cd ..

echo "Installing docs to web space..."

TOOLKIT_WEB_ROOT="/usr/global/web-pages/lc/www/toolkit"
DOCS_DIR_OLD="${TOOLKIT_WEB_ROOT}/docs_old"
DOCS_DIR="${TOOLKIT_WEB_ROOT}/docs"

echo "-----------------------------------------------------------------------"
if [ -d  ${DOCS_DIR_OLD} ]; then
    rm -rf ${DOCS_DIR_OLD}
    if [ $? -ne 0 ]; then
        echo "Error: 'rm' of docs_old failed"
        #exit 1
    fi
fi


if [ -d  ${DOCS_DIR} ]; then
    mv ${DOCS_DIR} ${DOCS_DIR_OLD}
    if [ $? -ne 0 ]; then
        echo "Error: 'mv' docs to docs_old failed"
        # exit 1
    fi
fi

cp -R ./${IP}/docs ${TOOLKIT_WEB_ROOT}/
if [ $? -ne 0 ]; then
    echo "Error: 'cp' failed"
    exit 1
fi

chgrp -R toolkit ${DOCS_DIR}
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
