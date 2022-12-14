#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/packages/cmake/cmake-3.16.8/bin/cmake
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: intel@19.0.4
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/usr/WS1/axom/libs/toss_3_x86_64_ib/2022_12_08_21_44_21/spack/lib/spack/env/intel/icc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/WS1/axom/libs/toss_3_x86_64_ib/2022_12_08_21_44_21/spack/lib/spack/env/intel/icpc" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/WS1/axom/libs/toss_3_x86_64_ib/2022_12_08_21_44_21/spack/lib/spack/env/intel/ifort" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/usr/tce/packages/intel/intel-19.0.4/bin/icc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/tce/packages/intel/intel-19.0.4/bin/icpc" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/tce/packages/intel/intel-19.0.4/bin/ifort" CACHE PATH "")

endif()

set(CMAKE_C_FLAGS "-gcc-name=/usr/tce/packages/gcc/gcc-8.1.0/bin/gcc" CACHE STRING "")

set(CMAKE_CXX_FLAGS "-gxx-name=/usr/tce/packages/gcc/gcc-8.1.0/bin/g++" CACHE STRING "")

set(CMAKE_Fortran_FLAGS "-gcc-name=/usr/tce/packages/gcc/gcc-8.1.0/bin/gcc" CACHE STRING "")

set(ENABLE_FORTRAN ON CACHE BOOL "")

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------

set(MPI_C_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3-intel-19.0.0/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3-intel-19.0.0/bin/mpicxx" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3-intel-19.0.0/bin/mpif90" CACHE PATH "")

set(MPIEXEC_EXECUTABLE "/usr/bin/srun" CACHE PATH "")

set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")

set(ENABLE_MPI ON CACHE BOOL "")

#------------------------------------------------------------------------------
# Hardware
#------------------------------------------------------------------------------

#------------------------------------------------
# Hardware Specifics
#------------------------------------------------

set(ENABLE_OPENMP OFF CACHE BOOL "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

#------------------------------------------------------------------------------
# TPLs
#------------------------------------------------------------------------------

# Root directory for generated TPLs

set(TPL_ROOT "/usr/WS1/axom/libs/toss_3_x86_64_ib/2022_12_08_21_44_21/intel-19.0.4" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.8.3" CACHE PATH "")

set(C2C_DIR "${TPL_ROOT}/c2c-1.3.0" CACHE PATH "")

set(MFEM_DIR "${TPL_ROOT}/mfem-4.5.0" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.22" CACHE PATH "")

set(LUA_DIR "/usr/WS1/axom/libs/toss_3_x86_64_ib/2022_12_08_21_44_21/clang-9.0.0/lua-5.4.4" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2022.03.0" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2022.03.1" CACHE PATH "")

set(CAMP_DIR "${TPL_ROOT}/camp-2022.03.2" CACHE PATH "")

# scr not built

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

# Root directory for generated developer tools

set(DEVTOOLS_ROOT "/collab/usr/gapps/axom/devtools/toss_3_x86_64_ib/2022_07_12_16_37_49/gcc-8.1.0" CACHE PATH "")

set(CLANGFORMAT_EXECUTABLE "/usr/tce/packages/clang/clang-10.0.0/bin/clang-format" CACHE PATH "")

set(PYTHON_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.9.13/bin/python3.9" CACHE PATH "")

set(ENABLE_DOCS ON CACHE BOOL "")

set(SPHINX_EXECUTABLE "${DEVTOOLS_ROOT}/py-sphinx-5.0.2/bin/sphinx-build" CACHE PATH "")

set(SHROUD_EXECUTABLE "${DEVTOOLS_ROOT}/py-shroud-0.12.2/bin/shroud" CACHE PATH "")

set(CPPCHECK_EXECUTABLE "${DEVTOOLS_ROOT}/cppcheck-2.8/bin/cppcheck" CACHE PATH "")

set(DOXYGEN_EXECUTABLE "${DEVTOOLS_ROOT}/doxygen-1.8.14/bin/doxygen" CACHE PATH "")


