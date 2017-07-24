################################
# documentation targets
################################

# add a target to generate documentation with Sphinx

add_sphinx_target( primal_docs )

if (DOXYGEN_FOUND)
  add_dependencies(primal_docs primal_doxygen)
endif()
