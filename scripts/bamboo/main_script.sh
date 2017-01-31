#!/bin/bash
# 09-12-2016 chang28, build-and-test.sh "clang@3.5.0" "Debug"
# 09-16-2016 chang28, build-and-test.sh "clang@3.5.0" "Debug" ""
# 09-19-2016 chang28, the decider has decided to have a configuration file call a main_script file, this is the main_script file, all environment variables are set up in the configuration file. 
# 01-30-2017 chang28, add -DENABLE_DOCS=false if $DOC= false

echo main_script version 0.9.5
echo "Configuring..."
echo "-----------------------------------------------------------------------"
if [ "$DOC" = false ]; then
   OPTIONS=$OPTIONS+" -DENABLE_DOCS=false";
fi
echo "Options: $OPTIONS"
./scripts/config-build.py $OPTIONS
if [ $? -ne 0 ]; then
    echo "Error: config-build.py failed"
    exit 1
fi
echo "-----------------------------------------------------------------------"

cd $BUILD_PATH

if [ "$BUILD" = true ]; then
    echo "Building..."
    echo "-----------------------------------------------------------------------"
    make VERBOSE=1 -j$JOBS
    if [ $? -ne 0 ]; then
        echo "Error: 'make' failed"
        exit 1
    fi
    echo "-----------------------------------------------------------------------"
fi

if [ "$TEST" = true ]; then
    echo "Running tests..."
    echo "-----------------------------------------------------------------------"
    make test ARGS="-T Test -j$JOBS"
    if [ $? -ne 0 ]; then
        echo "Error: 'make test' failed"
        exit 1
    fi
    echo "-----------------------------------------------------------------------"
fi

if [ "$DOC" = true ]; then
    echo "Making docs..."
    echo "-----------------------------------------------------------------------"
    make VERBOSE=1 docs
    if [ $? -ne 0 ]; then
        echo "Error: 'make docs' failed"
        exit 1
    fi
    echo "-----------------------------------------------------------------------"
fi

if [ "$INSTALL_FILES" = true ]; then
    echo "Installing files..."
    echo "-----------------------------------------------------------------------"
    make VERBOSE=1 install
    if [ $? -ne 0 ]; then
        echo "Error: 'make install' failed"
        exit 1
    fi
    echo "-----------------------------------------------------------------------"
fi

cd ..
if [ "$INSTALL_DOCS" = true ]; then


   echo "Installing docs to web space..."


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

   cp -R ./${INSTALL_PATH}/docs ${TOOLKIT_WEB_ROOT}/
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
fi

