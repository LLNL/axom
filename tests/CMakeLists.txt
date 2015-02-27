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