#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/bin/cmake
#------------------------------------------------------------------------------

set(CMAKE_PREFIX_PATH "/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_30_15_05_41/intel-2022.1.0/umpire-2023.06.0-6wmifjy5c2er4kkfqw7m5s4in7esiln5;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_30_15_05_41/intel-2022.1.0/raja-2023.06.0-urtfojwbjummyvyev7k4zn2oix53f4up;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_30_15_05_41/intel-2022.1.0/camp-2023.06.0-ovewmmz43udq34jrdlokh2zuzgkzrum4;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_30_15_05_41/intel-2022.1.0/mfem-4.5.2-dbvjel6jhwnwf6wzm76w6ghbjpy3v4gj;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_30_15_05_41/intel-2022.1.0/hypre-2.24.0-qdjyl5ibcpl4nxb6t5godfgvi5bb6hac;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_30_15_05_41/intel-2022.1.0/conduit-0.8.8-22refnw4twba4zznewcuczky2udimj7i;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_30_15_05_41/intel-2022.1.0/parmetis-4.0.3-nqxoj5gqogj32hfo5ognkcb6vg24zdlc;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_30_15_05_41/intel-2022.1.0/metis-5.1.0-tb5bijmooj2d6d2npjpajcj4fyu4yd4y;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_30_15_05_41/intel-2022.1.0/hdf5-1.8.22-6oeginq7kyyix35utjgsfxeghcnujtpg;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_30_15_05_41/intel-2022.1.0/c2c-1.8.0-gvlftgplgmtput3k4fy2nssecldmmlsa;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_30_15_05_41/intel-2022.1.0/gmake-4.4.1-e3r4ry6vgi7zp4mbwfokypkutqwpctek;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_30_15_05_41/intel-2022.1.0/blt-0.5.3-i7qltms7ksyz4xltznzbtsafywrykq7b;/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/latest/python-3.10.10;/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/latest/python-3.10.10;/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/latest/python-3.10.10;/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/latest/llvm-10.0.0;/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/latest/doxygen-1.9.6;/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/latest/cppcheck-2.9;/usr/tce/packages/mvapich2/mvapich2-2.3.6-intel-2022.1.0;/usr/tce" CACHE PATH "")

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "")

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: intel@=2022.1.0
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_30_15_05_41/spack/lib/spack/env/intel/icc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_30_15_05_41/spack/lib/spack/env/intel/icpc" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_30_15_05_41/spack/lib/spack/env/intel/ifort" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/usr/tce/packages/intel/intel-2022.1.0/compiler/2022.1.0/linux/bin/icx" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/tce/packages/intel/intel-2022.1.0/compiler/2022.1.0/linux/bin/icpx" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/tce/packages/intel/intel-2022.1.0/compiler/2022.1.0/linux/bin/ifx" CACHE PATH "")

endif()

set(CMAKE_GENERATOR "Unix Makefiles" CACHE STRING "")

set(ENABLE_FORTRAN ON CACHE BOOL "")

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------

set(MPI_C_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.6-intel-2022.1.0/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.6-intel-2022.1.0/bin/mpicxx" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.6-intel-2022.1.0/bin/mpif90" CACHE PATH "")

set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")

set(ENABLE_MPI ON CACHE BOOL "")

set(MPIEXEC_EXECUTABLE "/usr/bin/srun" CACHE PATH "")

#------------------------------------------------------------------------------
# Hardware
#------------------------------------------------------------------------------

#------------------------------------------------
# Hardware Specifics
#------------------------------------------------

set(ENABLE_OPENMP ON CACHE BOOL "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

#------------------------------------------------------------------------------
# TPLs
#------------------------------------------------------------------------------

set(TPL_ROOT "/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_30_15_05_41/intel-2022.1.0" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.8.8-22refnw4twba4zznewcuczky2udimj7i" CACHE PATH "")

set(C2C_DIR "${TPL_ROOT}/c2c-1.8.0-gvlftgplgmtput3k4fy2nssecldmmlsa" CACHE PATH "")

set(MFEM_DIR "${TPL_ROOT}/mfem-4.5.2-dbvjel6jhwnwf6wzm76w6ghbjpy3v4gj" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.22-6oeginq7kyyix35utjgsfxeghcnujtpg" CACHE PATH "")

set(LUA_DIR "/usr" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2023.06.0-urtfojwbjummyvyev7k4zn2oix53f4up" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2023.06.0-6wmifjy5c2er4kkfqw7m5s4in7esiln5" CACHE PATH "")

set(CAMP_DIR "${TPL_ROOT}/camp-2023.06.0-ovewmmz43udq34jrdlokh2zuzgkzrum4" CACHE PATH "")

# scr not built

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

set(DEVTOOLS_ROOT "/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/2023_10_17_16_15_32/._view/ypfidjpdm7fhkqcpekza67w5xgiaw7yg" CACHE PATH "")

set(CLANGFORMAT_EXECUTABLE "/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/latest/llvm-10.0.0/bin/clang-format" CACHE PATH "")

set(PYTHON_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/python3.10" CACHE PATH "")

set(ENABLE_DOCS ON CACHE BOOL "")

set(SPHINX_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/sphinx-build" CACHE PATH "")

set(SHROUD_EXECUTABLE "/collab/usr/gapps/shroud/public/toss_4_x86_64_ib/shroud-0.13.0/bin/shroud" CACHE PATH "")

set(CPPCHECK_EXECUTABLE "${DEVTOOLS_ROOT}/cppcheck-2.9/bin/cppcheck" CACHE PATH "")

set(DOXYGEN_EXECUTABLE "${DEVTOOLS_ROOT}/doxygen-1.9.6/bin/doxygen" CACHE PATH "")


