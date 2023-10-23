#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/bin/cmake
#------------------------------------------------------------------------------

set(CMAKE_PREFIX_PATH "/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_23_09_45_05/gcc-10.3.1/umpire-2023.06.0-z2ts74flz56i67azbbp5zyiwvwgwucrr;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_23_09_45_05/gcc-10.3.1/raja-2023.06.0-ca64lzglcihqwwjcmxuzfdvg7bhsd5rq;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_23_09_45_05/gcc-10.3.1/camp-2023.06.0-k4v6kc3zhnika36o2jvhncwiwcaifxqp;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_23_09_45_05/gcc-10.3.1/mfem-4.5.2-uyempb4k6nefxh36lxnbfb2dynpyaezr;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_23_09_45_05/gcc-10.3.1/hypre-2.24.0-lof7hdy7wsyuei436d6uhilprsgmr3ik;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_23_09_45_05/gcc-10.3.1/conduit-0.8.8-wngwwhpllk2gjewlizsk4plakvdu4beu;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_23_09_45_05/gcc-10.3.1/parmetis-4.0.3-lq6aryhur6qn3trc2qxizvz4nae5ngzf;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_23_09_45_05/gcc-10.3.1/metis-5.1.0-pyrzcrzqvl6q27kk637o2nwi3j7ibseb;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_23_09_45_05/gcc-10.3.1/hdf5-1.8.22-wikz2sxdo4zf76fdfl5do4mekkixf4yv;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_23_09_45_05/gcc-10.3.1/c2c-1.8.0-bsohtsh5a73tmmxlndgjuhzz6lzoif5j;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_23_09_45_05/gcc-10.3.1/gmake-4.4.1-y2kaxrqjbg3dcwo7dzfphztusfjqoowq;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_23_09_45_05/gcc-10.3.1/blt-0.5.3-iljv7wrfi4r5i4drs4yf65rgsuunaxjn;/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/latest/python-3.10.10;/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/latest/python-3.10.10;/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/latest/python-3.10.10;/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/latest/llvm-10.0.0;/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/latest/doxygen-1.9.6;/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/latest/cppcheck-2.9;/usr/tce/packages/mvapich2/mvapich2-2.3.6-gcc-10.3.1;/usr/tce;/usr/tce/packages/mvapich2/mvapich2-2.3.7-gcc-10.3.1" CACHE PATH "")

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "")

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: gcc@=10.3.1
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_23_09_45_05/spack/lib/spack/env/gcc/gcc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_23_09_45_05/spack/lib/spack/env/gcc/g++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_23_09_45_05/spack/lib/spack/env/gcc/gfortran" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/usr/tce/packages/gcc/gcc-10.3.1/bin/gcc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/tce/packages/gcc/gcc-10.3.1/bin/g++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/tce/packages/gcc/gcc-10.3.1/bin/gfortran" CACHE PATH "")

endif()

set(CMAKE_GENERATOR "Unix Makefiles" CACHE STRING "")

set(ENABLE_FORTRAN ON CACHE BOOL "")

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------

set(MPI_C_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.6-gcc-10.3.1/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.6-gcc-10.3.1/bin/mpicxx" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.6-gcc-10.3.1/bin/mpif90" CACHE PATH "")

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

set(TPL_ROOT "/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_23_09_45_05/gcc-10.3.1" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.8.8-wngwwhpllk2gjewlizsk4plakvdu4beu" CACHE PATH "")

set(C2C_DIR "${TPL_ROOT}/c2c-1.8.0-bsohtsh5a73tmmxlndgjuhzz6lzoif5j" CACHE PATH "")

set(MFEM_DIR "${TPL_ROOT}/mfem-4.5.2-uyempb4k6nefxh36lxnbfb2dynpyaezr" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.22-wikz2sxdo4zf76fdfl5do4mekkixf4yv" CACHE PATH "")

set(LUA_DIR "/usr" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2023.06.0-ca64lzglcihqwwjcmxuzfdvg7bhsd5rq" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2023.06.0-z2ts74flz56i67azbbp5zyiwvwgwucrr" CACHE PATH "")

set(CAMP_DIR "${TPL_ROOT}/camp-2023.06.0-k4v6kc3zhnika36o2jvhncwiwcaifxqp" CACHE PATH "")

# scr not built

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

set(DEVTOOLS_ROOT "/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/2023_10_17_16_15_32/._view/ypfidjpdm7fhkqcpekza67w5xgiaw7yg" CACHE PATH "")

set(CLANGFORMAT_EXECUTABLE "/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/latest/llvm-10.0.0/bin/clang-format" CACHE PATH "")

set(PYTHON_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/python3.10" CACHE PATH "")

set(ENABLE_DOCS ON CACHE BOOL "")

set(SPHINX_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/sphinx-build" CACHE PATH "")

set(SHROUD_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/shroud" CACHE PATH "")

set(CPPCHECK_EXECUTABLE "${DEVTOOLS_ROOT}/cppcheck-2.9/bin/cppcheck" CACHE PATH "")

set(DOXYGEN_EXECUTABLE "${DEVTOOLS_ROOT}/doxygen-1.9.6/bin/doxygen" CACHE PATH "")


