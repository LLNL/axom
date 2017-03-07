##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# toss_3_x86_64_ib-gcc@4.9.3
##################################

# cmake from uberenv
# cmake executable path: /usr/workspace/wsa/toolkit/thirdparty_libs/builds/2017_03_04_22_50_48/spack/opt/spack/toss_3_x86_64_ib/gcc-4.9.3/cmake-3.3.1-h62huv6j45n6zomdy76l23djjsai32yx/bin/cmake

#######
# using gcc@4.9.3 compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/usr/tce/packages/gcc/gcc-4.9.3/bin/gcc" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/tce/packages/gcc/gcc-4.9.3/bin/g++" CACHE PATH "")

# fortran compiler used by spack
set(ENABLE_FORTRAN ON CACHE PATH "")

set(CMAKE_Fortran_COMPILER  "/usr/tce/packages/gcc/gcc-4.9.3/bin/gfortran" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/2017_03_04_22_50_48/spack/opt/spack/toss_3_x86_64_ib/gcc-4.9.3/hdf5-1.8.16-xu3gjbihzhdtdoclaloujiac7xtuosy7" CACHE PATH "")

# conduit from uberenv
set(CONDUIT_DIR "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/2017_03_04_22_50_48/spack/opt/spack/toss_3_x86_64_ib/gcc-4.9.3/conduit-0.2.0-iezv674pvb6s3vzd24z577i3c5c74ocr" CACHE PATH "")

# doxygen from uberenv
set(DOXYGEN_EXECUTABLE "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/2017_03_04_22_50_48/spack/opt/spack/toss_3_x86_64_ib/gcc-4.9.3/doxygen-1.8.11-qvji5uglhzexrf2aocd6woqem5evpwlq/bin/doxygen" CACHE PATH "")

# python from uberenv
set(PYTHON_EXECUTABLE "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/2017_03_04_22_50_48/spack/opt/spack/toss_3_x86_64_ib/gcc-4.9.3/python-2.7.11-lea63zrc2ghlbkp23fwqihgaypyaunaa/bin/python" CACHE PATH "")

# lua from uberenv
set(LUA_DIR "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/2017_03_04_22_50_48/spack/opt/spack/toss_3_x86_64_ib/gcc-4.9.3/lua-5.1.5-esnntgauhb5eur2ftxv4cpxdrxq3wcd4" CACHE PATH "")

# sphinx from uberenv
set(SPHINX_EXECUTABLE "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/2017_03_04_22_50_48/spack/opt/spack/toss_3_x86_64_ib/gcc-4.9.3/python-2.7.11-lea63zrc2ghlbkp23fwqihgaypyaunaa/bin/sphinx-build" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/2017_03_04_22_50_48/spack/opt/spack/toss_3_x86_64_ib/gcc-4.9.3/uncrustify-0.61-gfugzvve5jxq4p6kpozhysf76qtcwxcc/bin/uncrustify" CACHE PATH "")

# sparsehash headers from uberenv
set(SPARSEHASH_DIR "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/2017_03_04_22_50_48/spack/opt/spack/toss_3_x86_64_ib/gcc-4.9.3/sparsehash-headers-2.0.2-czhsx5wtext3b4oipxtlk2bhmfaco7wr" CACHE PATH "")

# boost headers from uberenv
set(ENABLE_BOOST ON CACHE PATH "")
set(BOOST_ROOT "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/2017_03_04_22_50_48/spack/opt/spack/toss_3_x86_64_ib/gcc-4.9.3/boost-headers-1.58.0-nvmps5cyz6fnof343goaptduus4n65cz" CACHE PATH "")

# lcov and genhtml from uberenv
set(LCOV_PATH "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/2017_03_04_22_50_48/spack/opt/spack/toss_3_x86_64_ib/gcc-4.9.3/lcov-1.11-ms3k744sok6xi4ndnzcmwdxnhcqwpm3i/usr/bin/lcov" CACHE PATH "")

set(GENHTML_PATH "/usr/workspace/wsa/toolkit/thirdparty_libs/builds/2017_03_04_22_50_48/spack/opt/spack/toss_3_x86_64_ib/gcc-4.9.3/lcov-1.11-ms3k744sok6xi4ndnzcmwdxnhcqwpm3i/usr/bin/genhtml" CACHE PATH "")

##################################
# end uberenv host-config
##################################

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################
# Options added manually to 
# lc toss3 gcc@4.9.3  host configs
##############################################################################

##############################################################################
# MPI - manually added for now
##############################################################################
set("ENABLE_MPI" ON CACHE PATH "")
set("MPI_C_COMPILER" "/usr/tce/packages/mvapich2/mvapich2-2.2-gcc-4.9.3/bin/mpicc" CACHE PATH "")
set("MPI_CXX_COMPILER" "/usr/tce/packages/mvapich2/mvapich2-2.2-gcc-4.9.3/bin/mpiccxx" CACHE PATH "")
set("MPI_Fortran_COMPILER" "/usr/tce/packages/mvapich2/mvapich2-2.2-gcc-4.9.3/bin/mpif90" CACHE PATH "")

set("MPIEXEC" "/usr/bin/srun" CACHE PATH "")
set("MPIEXEC_NUMPROC_FLAG" "-n" CACHE PATH "")


##############################################################################
# !---------------------------------------------------------------------------
##############################################################################

