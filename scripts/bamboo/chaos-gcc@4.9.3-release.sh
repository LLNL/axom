#!/bin/sh
echo "Configuring..."
echo "-----------------------------------------------------------------------"
./scripts/config-build.py -c gcc@4.9.3 -buildtype Release --buildpath build-chaos-release
cd build-chaos-release
echo "-----------------------------------------------------------------------"

echo "Building..."
echo "-----------------------------------------------------------------------"
make -j16
echo "-----------------------------------------------------------------------"

echo "Running tests..."
echo "-----------------------------------------------------------------------"
make test ARGS="-T Test"
echo "-----------------------------------------------------------------------"

echo "Making docs..."
echo "-----------------------------------------------------------------------"
make docs
echo "-----------------------------------------------------------------------"

echo "Installing files..."
echo "-----------------------------------------------------------------------"
make install
echo "-----------------------------------------------------------------------"

echo "Installing docs to web space..."
echo "-----------------------------------------------------------------------"
#setenv ATK_JOB_BUILDDIR "/g/g16/atk/bambooAgent/asctoolkit.cab.llnl.gov/xml-data/build-dir/ASC-NIG-DC/asctoolkit/"
#/g/g16/atk/testing_framework/buildtest.py -s -b --builddir='/g/g16/atk/bambooAgent/asctoolkit.cab.llnl.gov/xml-data/build-dir/ASC-NIG-DC/asctoolkit/' --installdocs -t test docs install
/g/g16/atk/testing_framework/buildtest.py --installdocs --builddir=$ATK_JOB_BUILDDIR
cd  /usr/global/web-pages/lc/www/toolkit/
chgrp -R toolkit docs/
echo "----------------------
