#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/packages/cmake/cmake-3.21.1/bin/cmake
#------------------------------------------------------------------------------

set(CMAKE_PREFIX_PATH "/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/intel-19.0.4/umpire-2023.06.0-spxsyq3od3jpw7xarmaf27zkhpreokce;/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/intel-19.0.4/raja-2023.06.0-ddsf5gvw3o7fmxx6tqaboufen5rtq3z4;/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/intel-19.0.4/camp-2023.06.0-pl43tuxmbj3l3humx76qa2dnglwieb7u;/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/intel-19.0.4/mfem-4.5.2-gncr6owa7f72rqhpzooxegqm7mmpzrxo;/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/intel-19.0.4/hypre-2.24.0-3nrszj6z3zarr6ret4vogidku5l3dk4h;/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/intel-19.0.4/lua-5.4.4-zljt4isbv6foy6qelg5nbabt2qdrny22;/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/intel-19.0.4/ncurses-6.4-kdjhmjhd6mde3prcrx4n6yuon7ifsg6b;/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/intel-19.0.4/conduit-0.8.8-oe5nxzus2z6e4fkaig2qw6jjznvmbc54;/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/intel-19.0.4/parmetis-4.0.3-gll7fwyzerdprkjo6pfqtg53ucbtcwgh;/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/intel-19.0.4/metis-5.1.0-sbxpvwpfbuqmi7xahz577tegotyhxi7n;/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/intel-19.0.4/hdf5-1.8.22-q7z46mq5v45ezkpbskhjgzpovqwwjorl;/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/intel-19.0.4/c2c-1.8.0-a52qggtffivtdp7ef2zw5bdapzo3zg5o;/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/intel-19.0.4/gmake-4.4.1-mubwfgd2o6j6mfms6bkw5ppbzpitibqs;/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/intel-19.0.4/blt-0.5.3-ujtkyhrx64yk64oy2jpcsgxpbdz7tzsn;/collab/usr/gapps/axom/devtools/toss_3_x86_64_ib/latest/python-3.10.10;/collab/usr/gapps/axom/devtools/toss_3_x86_64_ib/latest/python-3.10.10;/collab/usr/gapps/axom/devtools/toss_3_x86_64_ib/latest/python-3.10.10;/usr/tce/packages/clang/clang-10.0.0;/collab/usr/gapps/axom/devtools/toss_3_x86_64_ib/latest/graphviz-7.1.0;/collab/usr/gapps/axom/devtools/toss_3_x86_64_ib/latest/doxygen-1.8.14;/collab/usr/gapps/axom/devtools/toss_3_x86_64_ib/latest/cppcheck-2.9;/usr/tce/packages/mvapich2/mvapich2-2.3-intel-19.0.0;/usr/tce/packages/cmake/cmake-3.21.1" CACHE PATH "")

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "")

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: intel@=19.0.4
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/spack/lib/spack/env/intel/icc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/spack/lib/spack/env/intel/icpc" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/spack/lib/spack/env/intel/ifort" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/usr/tce/packages/intel/intel-19.0.4/bin/icc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/tce/packages/intel/intel-19.0.4/bin/icpc" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/tce/packages/intel/intel-19.0.4/bin/ifort" CACHE PATH "")

endif()

set(CMAKE_C_FLAGS "-gcc-name=/usr/tce/packages/gcc/gcc-8.1.0/bin/gcc" CACHE STRING "")

set(CMAKE_CXX_FLAGS "-gxx-name=/usr/tce/packages/gcc/gcc-8.1.0/bin/g++" CACHE STRING "")

set(CMAKE_Fortran_FLAGS "-gcc-name=/usr/tce/packages/gcc/gcc-8.1.0/bin/gcc" CACHE STRING "")

set(CMAKE_GENERATOR "Unix Makefiles" CACHE STRING "")

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

set(TPL_ROOT "/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/intel-19.0.4" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.8.8-oe5nxzus2z6e4fkaig2qw6jjznvmbc54" CACHE PATH "")

set(C2C_DIR "${TPL_ROOT}/c2c-1.8.0-a52qggtffivtdp7ef2zw5bdapzo3zg5o" CACHE PATH "")

set(MFEM_DIR "${TPL_ROOT}/mfem-4.5.2-gncr6owa7f72rqhpzooxegqm7mmpzrxo" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.22-q7z46mq5v45ezkpbskhjgzpovqwwjorl" CACHE PATH "")

set(LUA_DIR "${TPL_ROOT}/lua-5.4.4-zljt4isbv6foy6qelg5nbabt2qdrny22" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2023.06.0-ddsf5gvw3o7fmxx6tqaboufen5rtq3z4" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2023.06.0-spxsyq3od3jpw7xarmaf27zkhpreokce" CACHE PATH "")

set(CAMP_DIR "${TPL_ROOT}/camp-2023.06.0-pl43tuxmbj3l3humx76qa2dnglwieb7u" CACHE PATH "")

# scr not built

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

set(DEVTOOLS_ROOT "/collab/usr/gapps/axom/devtools/toss_3_x86_64_ib/2023_04_18_13_40_46/._view/2axci4znbttcg4i676h7tlgjoffyysqt" CACHE PATH "")

set(CLANGFORMAT_EXECUTABLE "/usr/tce/packages/clang/clang-10.0.0/bin/clang-format" CACHE PATH "")

set(PYTHON_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/python3.10" CACHE PATH "")

set(ENABLE_DOCS ON CACHE BOOL "")

set(SPHINX_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/sphinx-build" CACHE PATH "")

set(SHROUD_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/shroud" CACHE PATH "")

set(CPPCHECK_EXECUTABLE "${DEVTOOLS_ROOT}/cppcheck-2.9/bin/cppcheck" CACHE PATH "")

set(DOXYGEN_EXECUTABLE "${DEVTOOLS_ROOT}/doxygen-1.8.14/bin/doxygen" CACHE PATH "")


