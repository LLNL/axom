blt_add_sphinx_target(lumberjack_docs)

if(DOXYGEN_FOUND)
   add_dependencies(lumberjack_docs lumberjack_doxygen) 
endif()
