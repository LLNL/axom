
################################
# documentation targets
################################

# add a target to generate documentation with Sphinx

if (BUILD_DOCS)
  add_sphinx_target( shroud_docs )
endif()
