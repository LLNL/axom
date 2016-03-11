##################################
# uberenv host-config
##################################
# chaos_5_x86_64_ib-gcc@4.7.1
##################################

# cmake from uberenv
# cmake exectuable path: /usr/gapps/asctoolkit/thirdparty_libs/2016_03_11/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.7.1/cmake-3.2.2-qoktsvcsmsw76d2hsnpgjjpnjwjovjov/bin/cmake

#######
# using gcc@4.7.1 compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/usr/apps/gnu/4.7.1/bin/gcc" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/apps/gnu/4.7.1/bin/g++" CACHE PATH "")

# fortran compiler used by spack
set(ENABLE_FORTRAN ON CACHE PATH "")

set(CMAKE_Fortran_COMPILER  "/usr/apps/gnu/4.7.1/bin/gfortran" CACHE PATH "")

# conduit from uberenv
set(CONDUIT_DIR "/usr/gapps/asctoolkit/thirdparty_libs/2016_03_11/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.7.1/conduit-github-mcb5wbbesfk75c3ixqzjajoycsw6cgah" CACHE PATH "")

# doxygen from uberenv
# doxygen exectuable path: /usr/gapps/asctoolkit/thirdparty_libs/2016_03_11/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.7.1/doxygen-1.8.10-tlyvcbyo2cigpqiji3r6ick43indqgip/bin/doxygen

# python from uberenv
set(PYTHON_EXECUTABLE "/usr/gapps/asctoolkit/thirdparty_libs/2016_03_11/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.7.1/python-2.7.8-g32izsb5gfir6ll5tsqa7hmcbvrpdblj/bin/python" CACHE PATH "")

# sphinx from uberenv
set(SPHINX_EXECUTABLE "/usr/gapps/asctoolkit/thirdparty_libs/2016_03_11/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.7.1/python-2.7.8-g32izsb5gfir6ll5tsqa7hmcbvrpdblj/bin/sphinx-build" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "/usr/gapps/asctoolkit/thirdparty_libs/2016_03_11/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.7.1/uncrustify-0.61-v65xp2l6i423faxv4oa7amaqykzjaa7e/bin/uncrustify" CACHE PATH "")

# sparsehash headers from uberenv
set(SPARSEHASH_DIR "/usr/gapps/asctoolkit/thirdparty_libs/2016_03_11/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.7.1/sparsehash-headers-2.0.2-ebexhuxsncs6fbkmof6djiv3d26rclik" CACHE PATH "")

# boost headers from uberenv
set(ENABLE_BOOST ON CACHE PATH "")
set(BOOST_ROOT "/usr/gapps/asctoolkit/thirdparty_libs/2016_03_11/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.7.1/boost-headers-1.58.0-zhkwg3db5a6xbucdytnjq3my52l5jlwu" CACHE PATH "")

# lcov and genhtml from uberenv
set(LCOV_PATH "/usr/gapps/asctoolkit/thirdparty_libs/2016_03_11/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.7.1/lcov-1.11-tvagnsw74h4kwgos22odgtzjcbegi2xt/usr/bin/lcov" CACHE PATH "")

set(GENHTML_PATH "/usr/gapps/asctoolkit/thirdparty_libs/2016_03_11/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.7.1/lcov-1.11-tvagnsw74h4kwgos22odgtzjcbegi2xt/usr/bin/genhtml" CACHE PATH "")

##################################
# end uberenv host-config
##################################

