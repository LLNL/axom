#-------------------------------------------------------------------------------
# Copyright (c) 2015, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# All rights reserved.
#
# This source code cannot be distributed without permission and further
# review from Lawrence Livermore National Laboratory.
#-------------------------------------------------------------------------------

set(quest_tests
    quest_boundingBox.cpp
    quest_bucket_tree.cpp
    quest_intersection.cpp
    quest_mesh_extent.cpp
    quest_morton.cpp
    quest_numeric_array.cpp
    quest_octree.cpp
    quest_orient.cpp
    quest_point.cpp
    quest_spatial_octree.cpp
    quest_squared_distance.cpp
    quest_vector.cpp
   )

foreach(test ${quest_tests})
  add_gtest( TEST_SOURCE ${test} DEPENDS_ON quest slic )
endforeach()





#------------------------------------------------------------------------------
#  The code below sets up a four tests for a clang bug reproducer
#  related to Morton indexing in 3D.
#
#  The tests are differentiated by different compiler defines
#  to ensure that each example is using exactly the same code.
#  These tests are setup to verify that that bug behavior is as expected.
#  The tests should fail if/when the bug is fixed in the clang compiler,
#  so we can remove the extra code related to the bugfix
#------------------------------------------------------------------------------

set(morton_reproducer_src    quest_morton_clang_repro.cpp)

# Tell the code if optimizations (-O2 and -O3 are enabled)
set(optsDefine "")
if(CMAKE_BUILD_TYPE MATCHES Release OR CMAKE_BUILD_TYPE MATCHES RelWithDebInfo)
  set(optsDefine "HAS_OPTS")
endif()

# Create a 2D version of the reproducer -- this should always produce the right results
make_executable( EXECUTABLE_NAME quest_mortonReproducer2D 
                 EXECUTABLE_SOURCE ${morton_reproducer_src} 
                 ADD_CTEST
               )
set_target_properties( quest_mortonReproducer2D
                        PROPERTIES COMPILE_DEFINITIONS "DIM=2;${optsDefine}" )

# Create the buggy 3D version of the reproducer
# It produces the wrong result, and the code tests for this
make_executable( EXECUTABLE_NAME quest_mortonReproducer3D_buggy1 
                 EXECUTABLE_SOURCE ${morton_reproducer_src}
                 ADD_CTEST 
               )
set_target_properties( quest_mortonReproducer3D_buggy1
                        PROPERTIES COMPILE_DEFINITIONS "DIM=3;${optsDefine}" )


# Create the still buggy 3D version of the reproducer -- with manual loop unrolling
# It produces the wrong result, and the code tests for this
make_executable( EXECUTABLE_NAME quest_mortonReproducer3D_buggy2 
                 EXECUTABLE_SOURCE ${morton_reproducer_src}
                 ADD_CTEST 
               )
set_target_properties( quest_mortonReproducer3D_buggy2
                        PROPERTIES COMPILE_DEFINITIONS "DIM=3;ATTEMPT_UNROLL;${optsDefine}" )

# Create the 'fixed'  3D version of the reproducer -- loop parameter is marked volatile
# The results should be correct this time.
make_executable( EXECUTABLE_NAME quest_mortonReproducer3D_bugFix
                 EXECUTABLE_SOURCE ${morton_reproducer_src} 
                 ADD_CTEST
               )
set_target_properties( quest_mortonReproducer3D_bugFix
                        PROPERTIES COMPILE_DEFINITIONS "DIM=3;APPLY_FIX;${optsDefine}" )


#------------------------------------------------------------------------------
