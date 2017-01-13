################################
# documentation targets
################################

add_subdirectory(web)

if(DOXYGEN_FOUND)
    add_subdirectory(doxygen)
endif()

if(SPHINX_FOUND)
    add_subdirectory(sphinx)
endif()
