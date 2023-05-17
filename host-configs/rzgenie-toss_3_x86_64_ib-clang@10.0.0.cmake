#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/packages/cmake/cmake-3.21.1/bin/cmake
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: clang@10.0.0
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_04_18_15_44_26/spack/lib/spack/env/clang/clang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_04_18_15_44_26/spack/lib/spack/env/clang/clang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_04_18_15_44_26/spack/lib/spack/env/clang/gfortran" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/usr/tce/packages/clang/clang-10.0.0/bin/clang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/tce/packages/clang/clang-10.0.0/bin/clang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/tce/packages/gcc/gcc-8.3.1/bin/gfortran" CACHE PATH "")

endif()

set(CMAKE_C_FLAGS "--gcc-toolchain=/usr/tce/packages/gcc/gcc-8.3.1" CACHE STRING "")

set(CMAKE_CXX_FLAGS "--gcc-toolchain=/usr/tce/packages/gcc/gcc-8.3.1" CACHE STRING "")

set(ENABLE_FORTRAN ON CACHE BOOL "")

set(BLT_EXE_LINKER_FLAGS " -Wl,-rpath,/usr/tce/packages/clang/clang-10.0.0/lib" CACHE STRING "Adds a missing libstdc++ rpath")

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------

set(MPI_C_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3-clang-10.0.0/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3-clang-10.0.0/bin/mpicxx" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3-clang-10.0.0/bin/mpif90" CACHE PATH "")

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

set(TPL_ROOT "/usr/WS1/axom/libs/toss_3_x86_64_ib/2023_04_18_15_44_26/clang-10.0.0" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.8.6-edn7xmhacopbvinkm4bhv2ebv24z45ix" CACHE PATH "")

set(C2C_DIR "${TPL_ROOT}/c2c-1.3.0-5lxac4njw5hmjtgm5i7wztlfnd4rqwiy" CACHE PATH "")

set(MFEM_DIR "${TPL_ROOT}/mfem-4.5.0-rozfjqgssgso4vg6fuaowokzrliyku5t" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.22-a2435y4n5jfwo2ujhxqtuto66anasm6r" CACHE PATH "")

set(LUA_DIR "${TPL_ROOT}/lua-5.4.4-6qbojxchmdentoxceahnanrzsm5wlfes" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2022.03.0-om5ba2dcyiomhp6fwfw6bfqw47v24aaa" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2022.03.1-s55j7me352jrw2qaytdqoqrwgw6ur56o" CACHE PATH "")

set(CAMP_DIR "${TPL_ROOT}/camp-2022.03.2-ptkto2bb4hfyaraqwpgklslyinukq4er" CACHE PATH "")

set(SCR_DIR "${TPL_ROOT}/scr-3.0.1-m4zltoufks4rcwsi6cyjf4fscvm3ndup" CACHE PATH "")

set(KVTREE_DIR "${TPL_ROOT}/kvtree-1.3.0-zgouicuhmkzoonhdbys3ysg45ydcnuan" CACHE PATH "")

set(DTCMP_DIR "${TPL_ROOT}/dtcmp-1.1.4-m5xpplpaobm2zmnyzeyvdlvsfp3uhyhm" CACHE PATH "")

set(SPATH_DIR "${TPL_ROOT}/spath-0.2.0-uzvf3jdoow32jphn52nr74f5yttpwubs" CACHE PATH "")

set(AXL_DIR "${TPL_ROOT}/axl-0.7.1-fd7aqekml3okahhkuc3lfdqt2x26heqy" CACHE PATH "")

set(LWGRP_DIR "${TPL_ROOT}/lwgrp-1.0.5-nzxc2ggdnkcf4fh6jqfgsevzg4ql3om3" CACHE PATH "")

set(ER_DIR "${TPL_ROOT}/er-0.2.0-5tzriyrrfxejvseuf7r4ubwyznpjwgku" CACHE PATH "")

set(RANKSTR_DIR "${TPL_ROOT}/rankstr-0.1.0-zghtcqxljvwoai4sivjv2amdlg36gr2o" CACHE PATH "")

set(REDSET_DIR "${TPL_ROOT}/redset-0.2.0-ufqtw7wpnzvx2dyutvre7tqyuiucy5fs" CACHE PATH "")

set(SHUFFILE_DIR "${TPL_ROOT}/shuffile-0.2.0-o5kwlrfg5ayslccnb6vdusnkrduexzkv" CACHE PATH "")

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


