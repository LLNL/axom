#
# Specify list of tests
#

set(catch_sidre_tests
   testCatchExample.cpp
   )


set(gtest_sidre_tests
#   testGTestExample.cpp
   sidre_smoke.cpp
   sidre_buffer.cpp
   sidre_view.cpp
   sidre_group.cpp
   sidre_opaque.cpp

   C_sidre_smoke.cpp
#   C_sidre_buffer.cpp
#   C_sidre_view.cpp
   C_sidre_group.cpp
#   C_sidre_opaque.cpp
   )



#
# Add catch based tests
#
foreach(test ${catch_sidre_tests})
    add_catch_test( TEST_SOURCE ${test} DEPENDS_ON sidre )
endforeach()

#
# Add gtest based tests
#
foreach(test ${gtest_sidre_tests})
    add_gtest( TEST_SOURCE ${test} DEPENDS_ON sidre conduit )
endforeach()


# sample fortran test
add_executable( testf testfortran.f )
set_target_properties(testf PROPERTIES Fortran_FORMAT "FREE")
target_link_libraries( testf sidre)
