#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/packages/cmake/cmake-3.21.1/bin/cmake
#------------------------------------------------------------------------------

set(CMAKE_PREFIX_PATH "/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/clang-9.0.0/umpire-2023.06.0-pxc3luewglsrjheqbmy4ethzzhtkp4wa;/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/clang-9.0.0/raja-2023.06.0-7s6l5urylci3jmqtttose3jarycrizpp;/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/clang-9.0.0/camp-2023.06.0-ojipssq33pvvq6qwllhxtb3rm7mygyby;/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/clang-9.0.0/mfem-4.5.2-os7vuxxmcxl2kwlxktmcyzfijg7nbcxn;/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/clang-9.0.0/hypre-2.24.0-k5mky734f76cd3wf4nsn7cu6ojb4qlng;/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/clang-9.0.0/lua-5.4.4-llfbtjxhew4bndeno5yiuhcsnhonmhx3;/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/clang-9.0.0/ncurses-6.4-vibranrcfja4nnd7363fxffvddp6cvan;/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/clang-9.0.0/conduit-0.8.8-adu2oxc6b3dgj2euult7r5jbqdhe4jsc;/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/clang-9.0.0/parmetis-4.0.3-feep6z2o6motj6ewaiyypzfixyvfuj6x;/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/clang-9.0.0/metis-5.1.0-gm5mheorgemz444kl4u452kpk37st3ip;/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/clang-9.0.0/hdf5-1.8.22-m6ch4naxrmq7eve3se5wkyetakirqp2x;/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/clang-9.0.0/c2c-1.8.0-l2mnatxn4h33nrbij5s5iq5fwy3egoch;/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/clang-9.0.0/gmake-4.4.1-7g2ykynv7s3hgjbm4qcxafa7ejhblyso;/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/clang-9.0.0/blt-0.5.3-au7nouiccwpfdujjkvk66k33g4d7fw3g;/collab/usr/gapps/axom/devtools/toss_3_x86_64_ib/latest/python-3.10.10;/collab/usr/gapps/axom/devtools/toss_3_x86_64_ib/latest/python-3.10.10;/collab/usr/gapps/axom/devtools/toss_3_x86_64_ib/latest/python-3.10.10;/usr/tce/packages/clang/clang-10.0.0;/collab/usr/gapps/axom/devtools/toss_3_x86_64_ib/latest/graphviz-7.1.0;/collab/usr/gapps/axom/devtools/toss_3_x86_64_ib/latest/doxygen-1.8.14;/collab/usr/gapps/axom/devtools/toss_3_x86_64_ib/latest/cppcheck-2.9;/usr/tce/packages/mvapich2/mvapich2-2.3-clang-9.0.0;/usr/tce/packages/cmake/cmake-3.21.1" CACHE PATH "")

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "")

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: clang@=9.0.0
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/spack/lib/spack/env/clang/clang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/spack/lib/spack/env/clang/clang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/spack/lib/spack/env/clang/gfortran" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/usr/tce/packages/clang/clang-9.0.0/bin/clang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/tce/packages/clang/clang-9.0.0/bin/clang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/tce/packages/gcc/gcc-8.3.1/bin/gfortran" CACHE PATH "")

endif()

set(CMAKE_C_FLAGS "--gcc-toolchain=/usr/tce/packages/gcc/gcc-8.3.1" CACHE STRING "")

set(CMAKE_CXX_FLAGS "--gcc-toolchain=/usr/tce/packages/gcc/gcc-8.3.1" CACHE STRING "")

set(CMAKE_GENERATOR "Unix Makefiles" CACHE STRING "")

set(ENABLE_FORTRAN ON CACHE BOOL "")

set(BLT_EXE_LINKER_FLAGS " -Wl,-rpath,/usr/tce/packages/clang/clang-9.0.0/lib" CACHE STRING "Adds a missing libstdc++ rpath")

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------

set(MPI_C_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3-clang-9.0.0/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3-clang-9.0.0/bin/mpicxx" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3-clang-9.0.0/bin/mpif90" CACHE PATH "")

set(MPIEXEC_EXECUTABLE "/usr/bin/srun" CACHE PATH "")

set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")

set(ENABLE_MPI ON CACHE BOOL "")

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

set(TPL_ROOT "/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_07_26_17_28_34/clang-9.0.0" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.8.8-adu2oxc6b3dgj2euult7r5jbqdhe4jsc" CACHE PATH "")

set(C2C_DIR "${TPL_ROOT}/c2c-1.8.0-l2mnatxn4h33nrbij5s5iq5fwy3egoch" CACHE PATH "")

set(MFEM_DIR "${TPL_ROOT}/mfem-4.5.2-os7vuxxmcxl2kwlxktmcyzfijg7nbcxn" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.22-m6ch4naxrmq7eve3se5wkyetakirqp2x" CACHE PATH "")

set(LUA_DIR "${TPL_ROOT}/lua-5.4.4-llfbtjxhew4bndeno5yiuhcsnhonmhx3" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2023.06.0-7s6l5urylci3jmqtttose3jarycrizpp" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2023.06.0-pxc3luewglsrjheqbmy4ethzzhtkp4wa" CACHE PATH "")

set(CAMP_DIR "${TPL_ROOT}/camp-2023.06.0-ojipssq33pvvq6qwllhxtb3rm7mygyby" CACHE PATH "")

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


