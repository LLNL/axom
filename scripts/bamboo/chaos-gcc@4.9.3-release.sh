#!/bin/sh
echo "Configuring..."
echo "-----------------------------------------------------------------------"
./scripts/config-build.py -c gcc@4.9.3 --buildtype Release --buildpath build-chaos-release --installpath install-chaos-release
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
# Is there a way to terminate if the make docs above failed?  So we don't try to copy over the
# files if that failed.
mv /usr/global/web-pages/lc/www/toolkit/docs/ /usr/global/web-pages/lc/www/toolkit/docs_old
rm -rf /usr/global/web-pages/lc/www/toolkit/docs
cp -R ./docs /usr/global/web-pages/lc/www/toolkit/
chgrp -R toolkit /usr/global/web-pages/lc/www/toolkit/docs/
echo "-----------------------------------------------------------------------"
