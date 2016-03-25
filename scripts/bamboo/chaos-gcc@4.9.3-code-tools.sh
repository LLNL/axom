#!/bin/sh
echo "Configuring..."
echo "-----------------------------------------------------------------------"
./scripts/config-build.py -c gcc@4.9.3 --buildtype Debug --cmakeoption="ENABLE_CODECOV=TRUE" 
cd build-chaos-gcc@4.9.3-debug
echo "-----------------------------------------------------------------------"

echo "-----------------------------------------------------------------------"
echo "Generating C/Fortran binding..."
make generate
echo "-----------------------------------------------------------------------"

echo "Building..."
echo "-----------------------------------------------------------------------"
make -j16
echo "-----------------------------------------------------------------------"

echo "Run tests and collect code coverage output."
echo "-----------------------------------------------------------------------"
make coverage ARGS="-T Test -E mpi"
echo "-----------------------------------------------------------------------"

echo "Install code coverage report to web space..."
echo "-----------------------------------------------------------------------"
pushd /usr/global/web-pages/lc/www/toolkit
rm -rf ./coverage
popd
cp -R ./coverage /usr/global/web-pages/lc/www/toolkit
echo "-----------------------------------------------------------------------"

echo "Running valgrind..."
echo "-----------------------------------------------------------------------"
make test ARGS=" -E mpi -D ExperimentalMemCheck --output-on-failure"
echo "-----------------------------------------------------------------------"

# TODO - Should we copy valgrind output to web area so team can access???
