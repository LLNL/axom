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
   sidre_class.cpp
   sidre_external.cpp
   sidre_group.cpp
   sidre_opaque.cpp
   sidre_view.cpp
   )

set(gtest_sidre_C_tests
   C_sidre_smoke.cpp
   C_sidre_buffer.cpp
   C_sidre_external.cpp
   C_sidre_group.cpp
   C_sidre_opaque.cpp
   C_sidre_types.cpp
   C_sidre_view.cpp
   )

set(fruit_sidre_tests
   F_sidre_smoke.f
   F_sidre_group.f
   F_sidre_buffer.f
   F_sidre_view.f
)



#
# Add catch based tests
#
foreach(test ${catch_sidre_tests})
    add_catch_test( TEST_SOURCE ${test} DEPENDS_ON sidre conduit slic)
endforeach()

#
# Add gtest C++ based tests
#
foreach(test ${gtest_sidre_tests})
    add_gtest( TEST_SOURCE ${test} DEPENDS_ON sidre conduit slic)
endforeach()

#
# Add gtest C based tests
#
if (ENABLE_FORTRAN)
  foreach(test ${gtest_sidre_C_tests})
      add_gtest( TEST_SOURCE ${test} DEPENDS_ON sidre conduit slic)
  endforeach()
endif()

#
# Add fortran tests
#
if (ENABLE_FORTRAN)
  foreach(test ${fruit_sidre_tests})
      add_fortran_test(
          TEST_SOURCE ${test}
          DEPENDS_ON fruit sidre slic conduit
      )
  endforeach()
endif()



