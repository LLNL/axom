################################
# MeshAPI
################################
project(meshapi)


################################
# Add the boost headers
################################
FIND_PACKAGE(Boost)
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
#add_subdirectory(examples)

################################
# Add tests
################################
if (BUILD_TESTING)
  add_subdirectory(tests)
endif()

################################
# Create CMake importable
# exports for all of our targets
################################
install(EXPORT ${PROJECT_NAME}-targets DESTINATION lib/cmake)
