##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# chaos_5_x86_64_ib-gcc@4.7.1
##################################

# cmake from uberenv
# cmake executable path: /usr/workspace/wsa/axom/thirdparty_libs/builds/2017_05_01_17_48_16/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.7.1/cmake-3.3.1-topto6ayaoqsg36mlt4tnxdafjtz3t2y/bin/cmake

#######
# using gcc@4.7.1 compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/usr/apps/gnu/4.7.1/bin/gcc" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/apps/gnu/4.7.1/bin/g++" CACHE PATH "")

# fortran compiler used by spack
set(ENABLE_FORTRAN "ON" CACHE PATH "")

set(CMAKE_Fortran_COMPILER "/usr/apps/gnu/4.7.1/bin/gfortran" CACHE PATH "")

# Root directory for generated TPLs
set(TPL_ROOT "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_05_01_17_48_16/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.7.1" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.16-xsesubbul5rxxvgxh2fo76qzxtn3qjgm" CACHE PATH "")

# conduit from uberenv
set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.2.1-3j42b5cmds7mv3nkrhra2m7uprclz5id" CACHE PATH "")

# boost headers from uberenv
set(BOOST_DIR "${TPL_ROOT}/boost-headers-1.58.0-zhkwg3db5a6xbucdytnjq3my52l5jlwu" CACHE PATH "")

# python from uberenv
set(PYTHON_EXECUTABLE "${TPL_ROOT}/python-2.7.11-mps32mjq56gsjox5ushhnt4ts5mio5sn/bin/python" CACHE PATH "")

# lua from uberenv
set(LUA_DIR "${TPL_ROOT}/lua-5.1.5-3m2omxvef6q3n6hfuimtcmajeeei2fa3" CACHE PATH "")

# doxygen from uberenv
set(DOXYGEN_EXECUTABLE "${TPL_ROOT}/doxygen-1.8.11-nzxxtekulhng4rxxjlji57wejhhkbxdm/bin/doxygen" CACHE PATH "")

# sphinx from uberenv
set(SPHINX_EXECUTABLE "${TPL_ROOT}/python-2.7.11-mps32mjq56gsjox5ushhnt4ts5mio5sn/bin/sphinx-build" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "${TPL_ROOT}/uncrustify-0.61-v65xp2l6i423faxv4oa7amaqykzjaa7e/bin/uncrustify" CACHE PATH "")

# lcov and genhtml from uberenv
set(LCOV_PATH "${TPL_ROOT}/lcov-1.11-tvagnsw74h4kwgos22odgtzjcbegi2xt/usr/bin/lcov" CACHE PATH "")

set(GENHTML_PATH "${TPL_ROOT}/lcov-1.11-tvagnsw74h4kwgos22odgtzjcbegi2xt/usr/bin/genhtml" CACHE PATH "")

##################################
# end uberenv host-config
##################################

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################
# Options added manually to 
# lc chaos5 gcc@4.7.1  host configs
##############################################################################

##############################################################################
# MPI - manually added for now.
##############################################################################
set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_HOME             "/usr/local/tools/mvapich2-gnu-2.0" CACHE PATH "")
set(MPI_C_COMPILER       "${MPI_HOME}/bin/mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER     "${MPI_HOME}/bin/mpicxx" CACHE PATH "")
set(MPI_Fortran_COMPILER "${MPI_HOME}/bin/mpif90" CACHE PATH "")

set(MPIEXEC "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

##############################################################################
# SHROUD - manually added for now. Use a public build add to TPL later
##############################################################################
set(SHROUD_EXECUTABLE "/usr/apps/shroud/bin/shroud" CACHE PATH "")

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################

