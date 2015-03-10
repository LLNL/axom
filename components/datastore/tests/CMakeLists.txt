#
# Specity list of tests
#

set(datastore_tests
   testExample.cxx
   )



#
# Build tests
#
foreach(test ${datastore_tests})
    add_catch_test( TEST_SOURCE ${test} DEPENDS_ON datastore )
endforeach()


## Add Google tests
add_gtest( TEST_SOURCE testSampleGTest.cxx DEPENDS_ON datastore )