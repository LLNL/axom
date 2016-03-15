##################################
# uberenv host-config
##################################
# chaos_5_x86_64_ib-clang@3.5.0
##################################

# cmake from uberenv
# cmake executable path: /usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/clang-3.5.0/cmake-3.3.1-2wgnfyfwvkzwxa7hhcsdju7t3jpx4fjp/bin/cmake

#######
# using clang@3.5.0 compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/usr/global/tools/clang/chaos_5_x86_64_ib/clang-omp-3.5.0/bin/clang" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/global/tools/clang/chaos_5_x86_64_ib/clang-omp-3.5.0/bin/clang++" CACHE PATH "")

# fortran compiler used by spack
# no fortran compiler

set(ENABLE_FORTRAN OFF CACHE PATH "")

# conduit from uberenv
set(CONDUIT_DIR "/usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/clang-3.5.0/conduit-github-rq4njsbhpfr4shjmxfqklzr5kbh4aiht" CACHE PATH "")

# doxygen from uberenv
set(DOXYGEN_EXECUTABLE "/usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/clang-3.5.0/doxygen-1.8.10-tmwreqbkrmnmaqnkxzumjanzmdevnn6r/bin/doxygen" CACHE PATH "")

# python from uberenv
set(PYTHON_EXECUTABLE "/usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/clang-3.5.0/python-2.7.8-lzn3frfwjy5akf5owasb74xiyq4qvhx6/bin/python" CACHE PATH "")

# sphinx from uberenv
set(SPHINX_EXECUTABLE "/usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/clang-3.5.0/python-2.7.8-lzn3frfwjy5akf5owasb74xiyq4qvhx6/bin/sphinx-build" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "/usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/clang-3.5.0/uncrustify-0.61-hfvykjq7rsj5v7kb3blntogdcuewppl2/bin/uncrustify" CACHE PATH "")

# sparsehash headers from uberenv
set(SPARSEHASH_DIR "/usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/clang-3.5.0/sparsehash-headers-2.0.2-afi2mzxg35xtfnrdfsci5toyjcwdolpq" CACHE PATH "")

# boost headers from uberenv
set(ENABLE_BOOST ON CACHE PATH "")
set(BOOST_ROOT "/usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/clang-3.5.0/boost-headers-1.58.0-a6n4lbbqdnhhkivpxazq5e4zrdhaspei" CACHE PATH "")

# lcov and genhtml from uberenv
set(LCOV_PATH "/usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/clang-3.5.0/lcov-1.11-otvm7uutci77rwj7cygg4d5vgecdorif/usr/bin/lcov" CACHE PATH "")

set(GENHTML_PATH "/usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/clang-3.5.0/lcov-1.11-otvm7uutci77rwj7cygg4d5vgecdorif/usr/bin/genhtml" CACHE PATH "")

##################################
# end uberenv host-config
##################################

