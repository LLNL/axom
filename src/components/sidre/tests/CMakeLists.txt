#
# Specify list of tests
#

set(catch_datastore_tests
   testCatchExample.cpp
   )


set(gtest_datastore_tests
   testGTestExample.cpp
   datastore_smoke.cpp
   datastore_buffer.cpp
   datastore_view.cpp
   datastore_group.cpp
   datastore_opaque.cpp
   )



#
# Add catch based tests
#
foreach(test ${catch_datastore_tests})
    add_catch_test( TEST_SOURCE ${test} DEPENDS_ON sidre )
endforeach()

#
# Add gtest based tests
#
foreach(test ${gtest_datastore_tests})
    add_gtest( TEST_SOURCE ${test} DEPENDS_ON sidre conduit )
endforeach()



