################################
# documentation targets
################################

# add a target to generate API documentation with Doxygen

add_custom_target(docs)

if(DOXYGEN_FOUND)
    add_subdirectory(doxygen)
endif()


if(SPHINX_FOUND)
#    add_subdirectory(sphinx)
endif()


