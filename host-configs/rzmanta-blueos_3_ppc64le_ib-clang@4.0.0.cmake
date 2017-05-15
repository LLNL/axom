##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# blueos_3_ppc64le_ib-clang@4.0.0
##################################

# cmake from uberenv
# cmake executable path: /usr/workspace/wsrzc/axom/thirdparty_libs/builds/2017_05_11_13_14_56/spack/opt/spack/blueos_3_ppc64le_ib/clang-4.0.0/cmake-3.3.1-pdq7ymiktz6yazhnfetwrnd7k2gqsby5/bin/cmake

#######
# using clang@4.0.0 compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/usr/tcetmp/packages/clang/clang-4.0.0/bin/clang" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/tcetmp/packages/clang/clang-4.0.0/bin/clang++" CACHE PATH "")

# fortran compiler used by spack
set(ENABLE_FORTRAN "ON" CACHE PATH "")

set(CMAKE_Fortran_COMPILER "/usr/tcetmp/packages/xl/xl-beta-2017.05.08/bin/xlf90" CACHE PATH "")

# Root directory for generated TPLs
set(TPL_ROOT "/usr/workspace/wsrzc/axom/thirdparty_libs/builds/2017_05_11_13_14_56/spack/opt/spack/blueos_3_ppc64le_ib/clang-4.0.0" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.16-oog4xh5irqbcv36y6rrxg46te5z4x7qm" CACHE PATH "")

# conduit from uberenv
set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.2.1-hmu3olzpd7ttg2fuzvxjgqwre7vmb6eb" CACHE PATH "")

# boost headers from uberenv
set(BOOST_DIR "${TPL_ROOT}/boost-headers-1.58.0-m3drhxrsm27ft72a6y55alsm2bcfmdgg" CACHE PATH "")

# python not build by uberenv

# lua not build by uberenv

# doxygen not built by uberenv

# sphinx not built by uberenv

# uncrustify not built by uberenv

# lcov and genhtml from uberenv
# lcov and genhtml not built by uberenv

##################################
# end uberenv host-config
##################################

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################
# Options added manually to 
# lc blueos clang@4.0.0  host configs
##############################################################################

##############################################################################
# MPI - manually added for now
##############################################################################

set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_HOME             "/usr/tcetmp/packages/spectrum-mpi/spectrum-mpi-2017.04.03-clang-coral-2017.05.03" CACHE PATH "")
set(MPI_Fortran_HOME         "/usr/tcetmp/packages/xl/xl-beta-2017.05.08" CACHE PATH "")
set(MPI_C_COMPILER           "${MPI_HOME}/bin/mpicc"   CACHE PATH "")
set(MPI_CXX_COMPILER         "${MPI_HOME}/bin/mpicxx"  CACHE PATH "")
set(MPI_Fortran_COMPILER     "${MPI_HOME}/bin/mpif77" CACHE PATH "")
set(MPI_Fortran_LIBRARIES    "${MPI_Fortran_HOME}/lib" CACHE PATH "")
set(MPI_Fortran_INCLUDE_PATH "${MPI_Fortran_HOME}/xlf/16.1.0/include" CACHE PATH "")

set(MPIEXEC              "mpirun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

##############################################################################
# SHROUD - manually added for now. Use a public build add to TPL later
##############################################################################
set(SHROUD_EXECUTABLE "/usr/apps/shroud/bin/shroud" CACHE PATH "")

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################

