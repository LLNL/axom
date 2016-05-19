#!/bin/bash

echo "Configuring..."
echo "-----------------------------------------------------------------------"
./scripts/config-build.py -c gcc@4.9.3 --buildtype Debug -DENABLE_CODECOV=TRUE 
if [ $? -ne 0 ]; then
    echo "Error: config-build.py failed"
fi
echo "-----------------------------------------------------------------------"

cd build-chaos-gcc@4.9.3-debug
    echo "Generating C/Fortran binding..."
    make generate
    if [ $? -ne 0 ]; then
        echo "Error: 'make generate' failed"
    fi
    echo "-----------------------------------------------------------------------"

    echo "Building..."
    echo "-----------------------------------------------------------------------"
    make -j16
    if [ $? -ne 0 ]; then
        echo "Error: 'make' failed"
    fi
    echo "-----------------------------------------------------------------------"

    echo "Run tests and collect code coverage output."
    echo "-----------------------------------------------------------------------"
    make coverage ARGS="-T Test -E mpi"
    if [ $? -ne 0 ]; then
        echo "Error: 'make coverage' failed"
    fi
    echo "-----------------------------------------------------------------------"

    echo "Install code coverage report to web space..."
    echo "-----------------------------------------------------------------------"
    pushd /usr/global/web-pages/lc/www/toolkit
        rm -rf ./coverage
        if [ $? -ne 0 ]; then
            echo "Error: 'rm' failed"
        fi
    popd
    cp -R ./coverage /usr/global/web-pages/lc/www/toolkit
    if [ $? -ne 0 ]; then
        echo "Error: 'cp' failed"
    fi
    echo "-----------------------------------------------------------------------"

    echo "Running valgrind..."
    echo "-----------------------------------------------------------------------"
    make test ARGS=" -E mpi -D ExperimentalMemCheck --output-on-failure"
    if [ $? -ne 0 ]; then
        echo "Error: 'make test' failed"
    fi
    echo "-----------------------------------------------------------------------"
cd ..
# TODO - Should we copy valgrind output to web area so team can access???
