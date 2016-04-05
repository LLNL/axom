# This script builds a minimal set of the toolkit components designed for integration testing of sidre in host application codes.
# This build configuration is intended to mimick how a host code may build the sidre component.

# This version is serial (no mpi).

echo "Configuring..."
echo "-----------------------------------------------------------------------"
./scripts/config-build.py -c gcc@4.9.3 --buildtype RelWithDebInfo -DBUILD_DOCS=OFF -DBUILD_EXAMPLES=OFF -DENABLE_quest=OFF -DENABLE_slam=OFF -DENABLE_shroud=OFF -DENABLE_PYTHON=OFF -DENABLE_MPI=OFF
cd build-chaos-gcc@4.9.3-relwithdebinfo
echo "-----------------------------------------------------------------------"

echo "Building..."
echo "-----------------------------------------------------------------------"
make -j16
echo "-----------------------------------------------------------------------"

echo "Run tests."
echo "-----------------------------------------------------------------------"
make test
echo "-----------------------------------------------------------------------"
