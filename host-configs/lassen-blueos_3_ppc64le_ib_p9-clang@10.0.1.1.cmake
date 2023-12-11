#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/packages/cmake/cmake-3.21.1/bin/cmake
#------------------------------------------------------------------------------

set(CMAKE_PREFIX_PATH "/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_18_12_46_33/clang-10.0.1.1/umpire-2023.06.0-h6wvcftrmg6w5e7etumawlsf5trgvsl7;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_18_12_46_33/clang-10.0.1.1/raja-2023.06.0-jra4z7ctml3mpzsuixtbegawg7g75cec;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_18_12_46_33/clang-10.0.1.1/camp-2023.06.0-abpqpyhwsuhzu4csmdy2nejqzrkfqh7a;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_18_12_46_33/clang-10.0.1.1/mfem-4.5.2-dbebdqverfh4lwgj5nal5d2vgmmvl4xb;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_18_12_46_33/clang-10.0.1.1/hypre-2.24.0-vrxf2vyvdlmylt5u2cu3vy4dqdc76yvq;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_18_12_46_33/clang-10.0.1.1/lua-5.4.4-o6kjgr3nt6ssvrj3mmonawxj3sbbcnel;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_18_12_46_33/clang-10.0.1.1/conduit-0.8.8-6vxxfa5atbct57263kwixtpa5jz54kew;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_18_12_46_33/clang-10.0.1.1/parmetis-4.0.3-aivfy5rovu5bwsiryry63apan5cdtc4m;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_18_12_46_33/clang-10.0.1.1/metis-5.1.0-k7t53hncr4qxzqoyzhekte5rqrorclyh;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_18_12_46_33/clang-10.0.1.1/hdf5-1.8.22-x626yhzvwnx7ssrmgl3lnvncx7s5nbkf;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_18_12_46_33/clang-10.0.1.1/zlib-1.2.13-aple3eiihvtvpurrdvosba4qab7rgdnr;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_18_12_46_33/clang-10.0.1.1/c2c-1.8.0-alkiietbi3stbihkzci37tecdzdyxpwh;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_18_12_46_33/clang-10.0.1.1/blt-0.5.3-e5rbowq4pvwwuchttguahzyjcpa5zvgr;/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/latest/python-3.10.10;/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/latest/python-3.10.10;/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/latest/python-3.10.10;/usr/tcetmp/packages/lapack/lapack-3.9.0-P9-gcc-7.3.1;/usr/tce/packages/clang/clang-10.0.0;/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/latest/graphviz-7.1.0;/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/latest/doxygen-1.9.6;/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/latest/cppcheck-2.9;/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-clang-10.0.1-gcc-8.3.1;/usr/tcetmp;/usr/tce/packages/cmake/cmake-3.21.1" CACHE PATH "")

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "")

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: clang@=10.0.1.1
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_18_12_46_33/spack/lib/spack/env/clang/clang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_18_12_46_33/spack/lib/spack/env/clang/clang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_18_12_46_33/spack/lib/spack/env/clang/gfortran" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/usr/tce/packages/clang/clang-ibm-10.0.1-gcc-8.3.1/bin/clang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/tce/packages/clang/clang-ibm-10.0.1-gcc-8.3.1/bin/clang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/tce/packages/gcc/gcc-8.3.1/bin/gfortran" CACHE PATH "")

endif()

set(CMAKE_C_STANDARD_LIBRARIES "-lgfortran" CACHE STRING "")

set(CMAKE_CXX_STANDARD_LIBRARIES "-lgfortran" CACHE STRING "")

set(CMAKE_Fortran_STANDARD_LIBRARIES "-lgfortran" CACHE STRING "")

set(CMAKE_GENERATOR "Unix Makefiles" CACHE STRING "")

set(ENABLE_FORTRAN ON CACHE BOOL "")

set(BLT_EXE_LINKER_FLAGS " -Wl,-rpath,/usr/tce/packages/clang/clang-ibm-10.0.1-gcc-8.3.1/lib" CACHE STRING "Adds a missing libstdc++ rpath")

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------

set(MPI_C_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-clang-10.0.1-gcc-8.3.1/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-clang-10.0.1-gcc-8.3.1/bin/mpicxx" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-clang-10.0.1-gcc-8.3.1/bin/mpif90" CACHE PATH "")

set(MPIEXEC_EXECUTABLE "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-clang-10.0.1-gcc-8.3.1/bin/mpirun" CACHE PATH "")

set(MPIEXEC_NUMPROC_FLAG "-np" CACHE STRING "")

set(ENABLE_MPI ON CACHE BOOL "")

set(BLT_MPI_COMMAND_APPEND "mpibind" CACHE STRING "")

#------------------------------------------------------------------------------
# Hardware
#------------------------------------------------------------------------------

#------------------------------------------------
# Hardware Specifics
#------------------------------------------------

set(ENABLE_OPENMP OFF CACHE BOOL "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

set(BLT_CMAKE_IMPLICIT_LINK_DIRECTORIES_EXCLUDE "/usr/tce/packages/gcc/gcc-4.9.3/lib64;/usr/tce/packages/gcc/gcc-4.9.3/lib64/gcc/powerpc64le-unknown-linux-gnu/4.9.3;/usr/tce/packages/gcc/gcc-4.9.3/gnu/lib64;/usr/tce/packages/gcc/gcc-4.9.3/gnu/lib64/gcc/powerpc64le-unknown-linux-gnu/4.9.3" CACHE STRING "")

#------------------------------------------------------------------------------
# TPLs
#------------------------------------------------------------------------------

set(TPL_ROOT "/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_18_12_46_33/clang-10.0.1.1" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.8.8-6vxxfa5atbct57263kwixtpa5jz54kew" CACHE PATH "")

set(C2C_DIR "${TPL_ROOT}/c2c-1.8.0-alkiietbi3stbihkzci37tecdzdyxpwh" CACHE PATH "")

set(MFEM_DIR "${TPL_ROOT}/mfem-4.5.2-dbebdqverfh4lwgj5nal5d2vgmmvl4xb" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.22-x626yhzvwnx7ssrmgl3lnvncx7s5nbkf" CACHE PATH "")

set(LUA_DIR "${TPL_ROOT}/lua-5.4.4-o6kjgr3nt6ssvrj3mmonawxj3sbbcnel" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2023.06.0-jra4z7ctml3mpzsuixtbegawg7g75cec" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2023.06.0-h6wvcftrmg6w5e7etumawlsf5trgvsl7" CACHE PATH "")

set(CAMP_DIR "${TPL_ROOT}/camp-2023.06.0-abpqpyhwsuhzu4csmdy2nejqzrkfqh7a" CACHE PATH "")

# scr not built

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

set(DEVTOOLS_ROOT "/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/2023_10_17_10_41_04/._view/5gug4et2qymmckk7yvtub4qcapp3figm" CACHE PATH "")

set(CLANGFORMAT_EXECUTABLE "/usr/tce/packages/clang/clang-10.0.0/bin/clang-format" CACHE PATH "")

set(PYTHON_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/python3.10" CACHE PATH "")

set(ENABLE_DOCS ON CACHE BOOL "")

set(SPHINX_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/sphinx-build" CACHE PATH "")

set(SHROUD_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/shroud" CACHE PATH "")

set(CPPCHECK_EXECUTABLE "${DEVTOOLS_ROOT}/cppcheck-2.9/bin/cppcheck" CACHE PATH "")

set(DOXYGEN_EXECUTABLE "${DEVTOOLS_ROOT}/doxygen-1.9.6/bin/doxygen" CACHE PATH "")


