################################
# documentation targets
################################

add_sphinx_target( mint_docs )

if ( DOXYGEN_FOUND )
  add_dependencies( mint_docs mint_doxygen )
endif()
