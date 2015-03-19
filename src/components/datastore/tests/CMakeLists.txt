#
# Specify list of tests
#

set(catch_datastore_tests
   testCatchExample.cxx
   )


set(gtest_datastore_tests
   testGTestExample.cxx
   datastore_smoke.cxx
   datastore_buffer.cxx
   )



#
# Add catch based tests
#
foreach(test ${catch_datastore_tests})
    add_catch_test( TEST_SOURCE ${test} DEPENDS_ON datastore )
endforeach()

#
# Add gtest based tests
#
foreach(test ${gtest_datastore_tests})
    add_gtest( TEST_SOURCE ${test} DEPENDS_ON datastore conduit)
endforeach()



