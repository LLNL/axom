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
   sidre_external.cpp
   sidre_group.cpp
   sidre_opaque.cpp
   sidre_types.cpp
   sidre_view.cpp
   C_sidre_smoke.cpp
   C_sidre_buffer.cpp
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


#
# Add fortran tests
#
add_fortran_test(TEST_SOURCE testfortran.f DEPENDS_ON sidre)




