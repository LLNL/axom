
################################
# documentation targets
################################

# add a target to generate documentation with Sphinx

if(DOXYGEN_FOUND)
  add_sphinx_target( sidre_docs )
  add_dependencies(sidre_docs sidre_doxygen) 
endif()

