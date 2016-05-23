#!/bin/bash

echo "Configuring..."
echo "-----------------------------------------------------------------------"
./scripts/config-build.py -c gcc@4.9.3 --buildtype Debug -DENABLE_COVERAGE=TRUE 
if [ $? -ne 0 ]; then
    echo "Error: config-build.py failed"
    exit 1
fi
echo "-----------------------------------------------------------------------"

cd build-chaos-gcc@4.9.3-debug
    echo "Generating C/Fortran binding..."
    make VERBOSE=1 generate
    if [ $? -ne 0 ]; then
        echo "Error: 'make generate' failed"
        exit 1
    fi
    echo "-----------------------------------------------------------------------"

    echo "Building..."
    echo "-----------------------------------------------------------------------"
    make VERBOSE=1 -j16
    if [ $? -ne 0 ]; then
        echo "Error: 'make' failed"
        exit 1
    fi
    echo "-----------------------------------------------------------------------"

    echo "Run tests and collect code coverage output."
    echo "-----------------------------------------------------------------------"
    make coverage ARGS="-T Test -E mpi"
    if [ $? -ne 0 ]; then
        echo "Error: 'make coverage' failed"
        exit 1
    fi
    echo "-----------------------------------------------------------------------"

    echo "Install code coverage report to web space..."
    echo "-----------------------------------------------------------------------"
    pushd /usr/global/web-pages/lc/www/toolkit
        rm -rf ./coverage
        if [ $? -ne 0 ]; then
            echo "Error: 'rm' failed"
            exit 1
        fi
    popd
    cp -R ./coverage /usr/global/web-pages/lc/www/toolkit
    if [ $? -ne 0 ]; then
        echo "Error: 'cp' failed"
        exit 1
    fi
    echo "-----------------------------------------------------------------------"

    echo "Running valgrind..."
    echo "-----------------------------------------------------------------------"
    make test ARGS=" -E mpi -D ExperimentalMemCheck --output-on-failure"
    if [ $? -ne 0 ]; then
        echo "Error: 'make test' failed"
        exit 1
    fi
    echo "-----------------------------------------------------------------------"
cd ..
# TODO - Should we copy valgrind output to web area so team can access???
