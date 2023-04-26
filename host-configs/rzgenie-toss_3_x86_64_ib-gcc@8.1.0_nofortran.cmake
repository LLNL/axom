#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/packages/cmake/cmake-3.21.1/bin/cmake
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: gcc@8.1.0
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_04_18_15_44_26/spack/lib/spack/env/gcc/gcc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_04_18_15_44_26/spack/lib/spack/env/gcc/g++" CACHE PATH "")

  # No Fortran compiler defined in spec
else()

  set(CMAKE_C_COMPILER "/usr/tce/packages/gcc/gcc-8.1.0/bin/gcc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/tce/packages/gcc/gcc-8.1.0/bin/g++" CACHE PATH "")

  # No Fortran compiler defined in spec
endif()

set(ENABLE_FORTRAN OFF CACHE BOOL "")

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------

set(MPI_C_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3-gcc-8.1.0/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3-gcc-8.1.0/bin/mpicxx" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3-gcc-8.1.0/bin/mpif90" CACHE PATH "")

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

# Root directory for generated TPLs

set(TPL_ROOT "/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_04_18_15_44_26/gcc-8.1.0" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.8.6-tfjlaocjusulywyxjkosrmeekj27sock" CACHE PATH "")

set(C2C_DIR "${TPL_ROOT}/c2c-1.3.0-5hswnx652pz3mp3izgtw7nvjxhhlu5c3" CACHE PATH "")

# MFEM not built

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.22-4mdspjh2lyclwkeyjzfaoduv4sknfldu" CACHE PATH "")

set(LUA_DIR "${TPL_ROOT}/lua-5.4.4-44id5e5sxu3n6k2ikn3xuk3wx44db4vk" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2022.03.0-ragtq6i2m3nkzj5gubrikmeg3ummqebs" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2022.03.1-vmcr7qr2bhqfrmhdw3r7yctkcqusqqo7" CACHE PATH "")

set(CAMP_DIR "${TPL_ROOT}/camp-2022.03.2-ht2ydkfs4pznfxyzeqh2wm26ucdvla6w" CACHE PATH "")

set(SCR_DIR "${TPL_ROOT}/scr-3.0.1-w4u3tgvwhpxoa3ckvjvgcctvhuln4b6v" CACHE PATH "")

set(KVTREE_DIR "${TPL_ROOT}/kvtree-1.3.0-sawpcgcgvg6momllj7iy5xb2ro23mexn" CACHE PATH "")

set(DTCMP_DIR "${TPL_ROOT}/dtcmp-1.1.4-cjeq3puv72cnoivogbkivtivyo6ml3u2" CACHE PATH "")

set(SPATH_DIR "${TPL_ROOT}/spath-0.2.0-nfdrf7qdkofrjkspcahw2lidz7ylgj6w" CACHE PATH "")

set(AXL_DIR "${TPL_ROOT}/axl-0.7.1-bchhoq56i2amef3w7y6lj6n7polalki5" CACHE PATH "")

set(LWGRP_DIR "${TPL_ROOT}/lwgrp-1.0.5-yvycp3ob6io3xyz7h53omx4br5mduh76" CACHE PATH "")

set(ER_DIR "${TPL_ROOT}/er-0.2.0-rjl2yhe6pp7orn3u3z7ntqcmqhzjscpt" CACHE PATH "")

set(RANKSTR_DIR "${TPL_ROOT}/rankstr-0.1.0-npxlcaizrs5p7pqykosfljxn2mole5sp" CACHE PATH "")

set(REDSET_DIR "${TPL_ROOT}/redset-0.2.0-etwtxo5bm4bhjaw2k7i3pqlds4aoiavm" CACHE PATH "")

set(SHUFFILE_DIR "${TPL_ROOT}/shuffile-0.2.0-5vinuidfvzwpwvlixzjfld34fvfeifgy" CACHE PATH "")

set(LIBYOGRT_DIR "/usr" CACHE PATH "")

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

# Root directory for generated developer tools

set(DEVTOOLS_ROOT "/collab/usr/gapps/axom/devtools/toss_3_x86_64_ib/2023_04_18_13_40_46/._view/2axci4znbttcg4i676h7tlgjoffyysqt" CACHE PATH "")

set(CLANGFORMAT_EXECUTABLE "/usr/tce/packages/clang/clang-10.0.0/bin/clang-format" CACHE PATH "")

set(PYTHON_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/python3.10" CACHE PATH "")

set(ENABLE_DOCS ON CACHE BOOL "")

set(SPHINX_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/sphinx-build" CACHE PATH "")

set(SHROUD_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/shroud" CACHE PATH "")

set(CPPCHECK_EXECUTABLE "${DEVTOOLS_ROOT}/cppcheck-2.9/bin/cppcheck" CACHE PATH "")

set(DOXYGEN_EXECUTABLE "${DEVTOOLS_ROOT}/doxygen-1.8.14/bin/doxygen" CACHE PATH "")


