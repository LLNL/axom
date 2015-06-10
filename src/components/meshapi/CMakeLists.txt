################################
# MeshAPI
################################
project(meshapi)



#SET(CXX_11 "")     
#SET(CXX_11 "-std=c++0x -Wno-deprecated-declarations")     

### -Wno-missing-braces added to avoid warnings from boost
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_11} -Wall -Wno-missing-braces -fno-strict-aliasing " CACHE INTERNAL "meshapi") # -Wextra

### Need to add some extra flags to allow parameters and variables that are only used in ATK_WARNINGs
SET(EXTRA_FLAGS_WHEN_NO_DEBUG "-Wno-unused-parameter -Wno-unused-variable")
SET(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} ${EXTRA_FLAGS_WHEN_NO_DEBUG}" CACHE INTERNAL "meshapi") # -Wextra
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${EXTRA_FLAGS_WHEN_NO_DEBUG}" CACHE INTERNAL "meshapi") # -Wextra
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${EXTRA_FLAGS_WHEN_NO_DEBUG}" CACHE INTERNAL "meshapi") # -Wextra


#SET(CMAKE_CXX_VECTORIZER_FLAGS "-ftree-vectorize -fdump-tree-vect -ftree-vectorizer-verbose=7 -fdump-tree-optimized-details ") #-march=native "

#SET(CMAKE_CXX_FLAGS_DEBUG "-g -O0 ")
#SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -g ${CMAKE_CXX_VECTORIZER_FLAGS}")  #"-g0 -O3 -v"


################################
# Add the boost headers
################################
#FIND_PACKAGE(Boost)
MESSAGE(STATUS "Boost include dir: " ${Boost_INCLUDE_DIR})
MESSAGE(STATUS "Boost libraries dir: " ${BOOST_LIBRARY_DIR} )
MESSAGE(STATUS "Boost version: " ${BOOST_VERSION} )

################################
# Add the mesh api sources
################################
add_subdirectory(src)

################################
# Add examples
################################
add_subdirectory(examples)

################################
# Add tests
################################
if (BUILD_TESTING)
  add_subdirectory(tests)
endif()

add_code_check_targets(uncrustify.cfg)

################################
# Add docs
################################
add_subdirectory(docs)


################################
# Create CMake importable
# exports for all of our targets
################################
install(EXPORT ${PROJECT_NAME}-targets DESTINATION lib/cmake)
