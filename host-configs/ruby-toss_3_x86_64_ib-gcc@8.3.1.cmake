#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/packages/cmake/cmake-3.21.1/bin/cmake
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: gcc@8.3.1
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_04_18_17_37_32/spack/lib/spack/env/gcc/gcc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_04_18_17_37_32/spack/lib/spack/env/gcc/g++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_04_18_17_37_32/spack/lib/spack/env/gcc/gfortran" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/usr/tce/packages/gcc/gcc-8.3.1/bin/gcc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/tce/packages/gcc/gcc-8.3.1/bin/g++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/tce/packages/gcc/gcc-8.3.1/bin/gfortran" CACHE PATH "")

endif()

set(ENABLE_FORTRAN ON CACHE BOOL "")

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------

set(MPI_C_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3-gcc-8.3.1/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3-gcc-8.3.1/bin/mpicxx" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3-gcc-8.3.1/bin/mpif90" CACHE PATH "")

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

set(TPL_ROOT "/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_04_18_17_37_32/gcc-8.3.1" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.8.6-jjokus6cknvjvqzdzvt64yqij2omfopx" CACHE PATH "")

set(C2C_DIR "${TPL_ROOT}/c2c-1.3.0-njdbin3s4wxyoq3a7i7nphp6gc4pzong" CACHE PATH "")

set(MFEM_DIR "${TPL_ROOT}/mfem-4.5.0-oatd6hooasqxqyxmkrp3xtdjtik6omdn" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.22-jjcf433zhyg7nyius6krvmqh6gtubozs" CACHE PATH "")

set(LUA_DIR "${TPL_ROOT}/lua-5.4.4-fohlqsc25bnexnmtysardxoh3ufob3ek" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2022.03.0-nh56ughqsq4muuqeetuqxnisdow3jajt" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2022.03.1-fqwvpqlnu5tqeonpsvkko6gn2kad4lsc" CACHE PATH "")

set(CAMP_DIR "${TPL_ROOT}/camp-2022.03.2-t66aoj353nnv3vejfvxemzcicprckaxy" CACHE PATH "")

set(SCR_DIR "${TPL_ROOT}/scr-3.0.1-7lhdj2k5gtv74so2vxk7qyijg6wwsfxc" CACHE PATH "")

set(KVTREE_DIR "${TPL_ROOT}/kvtree-1.3.0-oh5uxw354uone2djvxuenai2fkt7wgk2" CACHE PATH "")

set(DTCMP_DIR "${TPL_ROOT}/dtcmp-1.1.4-4pfp5n4xxtkdx243ln5yvbl3mrgymchd" CACHE PATH "")

set(SPATH_DIR "${TPL_ROOT}/spath-0.2.0-nivjuuwu67ae5dd2lxlespntemurpbsj" CACHE PATH "")

set(AXL_DIR "${TPL_ROOT}/axl-0.7.1-eu4edsshjqckeeiqw5jm3mv4s7clyrhc" CACHE PATH "")

set(LWGRP_DIR "${TPL_ROOT}/lwgrp-1.0.5-i4okl7zfpdz4nt6vmqbn5hdqu7sj4zix" CACHE PATH "")

set(ER_DIR "${TPL_ROOT}/er-0.2.0-fovnjpn4ohrgrcje7lmeordgl673ipcq" CACHE PATH "")

set(RANKSTR_DIR "${TPL_ROOT}/rankstr-0.1.0-vfzhd3wjqqczidapi6c7oj2qcyr7bxib" CACHE PATH "")

set(REDSET_DIR "${TPL_ROOT}/redset-0.2.0-ywzpb7q6leswuotrmkwty5qg4becfdk5" CACHE PATH "")

set(SHUFFILE_DIR "${TPL_ROOT}/shuffile-0.2.0-2urru52aqe2ecn6bf2o57zdgaqann6zo" CACHE PATH "")

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


