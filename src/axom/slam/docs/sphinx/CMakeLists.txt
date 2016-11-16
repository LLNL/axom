
################################
# documentation targets
################################

# add a target to generate documentation with Sphinx

add_sphinx_target( slam_docs )

if(DOXYGEN_FOUND)
   add_dependencies(slam_docs slam_doxygen) 
endif()

