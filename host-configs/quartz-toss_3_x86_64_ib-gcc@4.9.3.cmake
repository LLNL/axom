##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# toss_3_x86_64_ib-gcc@4.9.3
##################################

# cmake from uberenv
# cmake executable path: /usr/workspace/wsa/axom/thirdparty_libs/builds/2017_05_01_17_47_34/spack/opt/spack/toss_3_x86_64_ib/gcc-4.9.3/cmake-3.3.1-h62huv6j45n6zomdy76l23djjsai32yx/bin/cmake

#######
# using gcc@4.9.3 compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/usr/tce/packages/gcc/gcc-4.9.3/bin/gcc" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/tce/packages/gcc/gcc-4.9.3/bin/g++" CACHE PATH "")

# fortran compiler used by spack
set(ENABLE_FORTRAN "ON" CACHE PATH "")

set(CMAKE_Fortran_COMPILER "/usr/tce/packages/gcc/gcc-4.9.3/bin/gfortran" CACHE PATH "")

# Root directory for generated TPLs
set(TPL_ROOT "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_05_01_17_47_34/spack/opt/spack/toss_3_x86_64_ib/gcc-4.9.3" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.16-v4jrbjlu2s6vl5a67kv3ydg6ssdsi2tj" CACHE PATH "")

# conduit from uberenv
set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.2.1-aucg4njv6djcvymsftkrzwswblstl7hs" CACHE PATH "")

# boost headers from uberenv
set(BOOST_DIR "${TPL_ROOT}/boost-headers-1.58.0-nvmps5cyz6fnof343goaptduus4n65cz" CACHE PATH "")

# python from uberenv
set(PYTHON_EXECUTABLE "${TPL_ROOT}/python-2.7.11-lea63zrc2ghlbkp23fwqihgaypyaunaa/bin/python" CACHE PATH "")

# lua from uberenv
set(LUA_DIR "${TPL_ROOT}/lua-5.1.5-esnntgauhb5eur2ftxv4cpxdrxq3wcd4" CACHE PATH "")

# doxygen from uberenv
set(DOXYGEN_EXECUTABLE "${TPL_ROOT}/doxygen-1.8.11-qvji5uglhzexrf2aocd6woqem5evpwlq/bin/doxygen" CACHE PATH "")

# sphinx from uberenv
set(SPHINX_EXECUTABLE "${TPL_ROOT}/python-2.7.11-lea63zrc2ghlbkp23fwqihgaypyaunaa/bin/sphinx-build" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "${TPL_ROOT}/uncrustify-0.61-gfugzvve5jxq4p6kpozhysf76qtcwxcc/bin/uncrustify" CACHE PATH "")

# lcov and genhtml from uberenv
set(LCOV_PATH "${TPL_ROOT}/lcov-1.11-ms3k744sok6xi4ndnzcmwdxnhcqwpm3i/usr/bin/lcov" CACHE PATH "")

set(GENHTML_PATH "${TPL_ROOT}/lcov-1.11-ms3k744sok6xi4ndnzcmwdxnhcqwpm3i/usr/bin/genhtml" CACHE PATH "")

##################################
# end uberenv host-config
##################################

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################
# Options added manually to 
# lc toss3 gcc@4.9.3  host configs
##############################################################################

##############################################################################
# MPI - manually added for now
##############################################################################
set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_HOME             "/usr/tce/packages/mvapich2/mvapich2-2.2-gcc-4.9.3" CACHE PATH "")
set(MPI_C_COMPILER       "${MPI_HOME}/bin/mpicc"   CACHE PATH "")
set(MPI_CXX_COMPILER     "${MPI_HOME}/bin/mpicxx"  CACHE PATH "")
set(MPI_Fortran_COMPILER "${MPI_HOME}/bin/mpifort" CACHE PATH "")

set(MPIEXEC              "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

##############################################################################
# SHROUD - manually added for now. Use a public build add to TPL later
##############################################################################
set(SHROUD_EXECUTABLE "/usr/apps/shroud/bin/shroud" CACHE PATH "")

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################

