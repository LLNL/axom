
################################
# documentation targets
################################

# add a target to generate documentation with Sphinx

add_sphinx_target( datastore_docs )

if(DOXYGEN_FOUND)
   add_dependencies(datastore_docs datastore_doxygen) 
endif()

