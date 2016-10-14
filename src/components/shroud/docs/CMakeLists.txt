
################################
# documentation targets
################################

# add a target to generate documentation with Sphinx
if(SPHINX_FOUND)
    add_sphinx_target( shroud_docs )
endif()