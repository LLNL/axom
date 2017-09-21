
################################
# documentation targets
################################

# add a target to generate documentation with Sphinx

blt_add_sphinx_target( spio_docs )

if(DOXYGEN_FOUND)
   add_dependencies(spio_docs spio_doxygen)
endif()


