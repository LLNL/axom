#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/bin/cmake
#------------------------------------------------------------------------------

set(BLT_CXX_STD "c++17" CACHE STRING "")

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "")

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: intel@=2022.1.0
#------------------------------------------------------------------------------

set(CMAKE_C_COMPILER "/usr/tce/packages/intel/intel-2022.1.0-magic/bin/icx" CACHE PATH "")

set(CMAKE_CXX_COMPILER "/usr/tce/packages/intel/intel-2022.1.0-magic/bin/icpx" CACHE PATH "")

set(CMAKE_Fortran_COMPILER "/usr/tce/packages/intel/intel-2022.1.0-magic/bin/ifx" CACHE PATH "")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --gcc-toolchain=/usr/tce/packages/gcc/gcc-12.1.1-magic" CACHE STRING "")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --gcc-toolchain=/usr/tce/packages/gcc/gcc-12.1.1-magic" CACHE STRING "")

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} --gcc-toolchain=/usr/tce/packages/gcc/gcc-12.1.1-magic" CACHE STRING "")

set(CMAKE_GENERATOR "Unix Makefiles" CACHE STRING "")

set(ENABLE_FORTRAN OFF CACHE BOOL "")

set(ENABLE_EXAMPLES OFF CACHE BOOL "")
#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------

set(MPI_C_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.7-intel-2022.1.0-magic/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.7-intel-2022.1.0-magic/bin/mpicxx" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.7-intel-2022.1.0-magic/bin/mpif90" CACHE PATH "")

set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")

set(ENABLE_MPI ON CACHE BOOL "")

set(MPIEXEC_EXECUTABLE "/usr/bin/srun" CACHE PATH "")

#------------------------------------------------------------------------------
# Hardware
#------------------------------------------------------------------------------

#------------------------------------------------
# Hardware Specifics
#------------------------------------------------

set(ENABLE_OPENMP OFF CACHE BOOL "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

set(AXOM_ENABLE_INLET OFF CACHE BOOL "")
set(AXOM_ENABLE_KLEE OFF CACHE BOOL "")
set(AXOM_ENABLE_SIDRE OFF CACHE BOOL "")
set(AXOM_ENABLE_TOOLS OFF CACHE BOOL "")

#------------------------------------------------------------------------------
# TPLs
#------------------------------------------------------------------------------

set(CONDUIT_DIR "" CACHE PATH "")

set(C2C_DIR "" CACHE PATH "")
set(AXOM_USE_C2C OFF CACHE BOOL "")

set(MFEM_DIR "" CACHE PATH "")
set(AXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION OFF CACHE BOOL "")

set(HDF5_DIR "/usr/WS1/ale3d_au/gapps/ale3d/toss_4_x86_64_ib/intel/PUBLIC/packages/hdf5/1.8.10.1f/compiler-intel-2022.1.0-gcc-12.1.1/cpp-17/cmake-3.21.1/zlib-1.2.11.b/share/cmake/hdf5" CACHE PATH "")

set(LUA_DIR "/usr" CACHE PATH "")

set(RAJA_DIR "/usr/WS1/ale3d_au/gapps/ale3d/toss_4_x86_64_ib/intel/PUBLIC/packages/raja/2022.10.5.c/compiler-intel-2022.1.0-gcc-12.1.1/cpp-17/blt-0.5.2.c/cmake-3.21.1/camp-2022.10.1.b/cuda-none/lib/cmake/raja" CACHE PATH "")

set(UMPIRE_DIR "/usr/WS1/ale3d_au/gapps/ale3d/toss_4_x86_64_ib/intel/PUBLIC/packages/umpire/2022.10.0.a/compiler-intel-2022.1.0-gcc-12.1.1/cpp-17/blt-0.5.2.c/cmake-3.21.1/camp-2022.10.1.b/cuda-none/lib/cmake/umpire" CACHE PATH "")

set(CAMP_DIR "/usr/WS1/ale3d_au/gapps/ale3d/toss_4_x86_64_ib/intel/PUBLIC/packages/camp/2022.10.1.b/compiler-intel-2022.1.0-gcc-12.1.1/cpp-17/blt-0.5.2.c/cmake-3.21.1/cuda-none/lib/cmake/camp" CACHE PATH "")

# scr not built

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

if (0)
set(DEVTOOLS_ROOT "/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/2023_10_17_16_15_32/._view/ypfidjpdm7fhkqcpekza67w5xgiaw7yg" CACHE PATH "")

set(CLANGFORMAT_EXECUTABLE "/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/latest/llvm-10.0.0/bin/clang-format" CACHE PATH "")

set(PYTHON_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/python3.10" CACHE PATH "")

set(ENABLE_DOCS OFF CACHE BOOL "")

set(SPHINX_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/sphinx-build" CACHE PATH "")

set(SHROUD_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/shroud" CACHE PATH "")

set(CPPCHECK_EXECUTABLE "${DEVTOOLS_ROOT}/cppcheck-2.9/bin/cppcheck" CACHE PATH "")

set(DOXYGEN_EXECUTABLE "${DEVTOOLS_ROOT}/doxygen-1.9.6/bin/doxygen" CACHE PATH "")
endif()
