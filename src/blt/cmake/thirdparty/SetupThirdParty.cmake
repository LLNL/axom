####################################
# BLT 3rd Party Lib Support
####################################

################################
# Documentation Packages
################################
if (DOXYGEN_EXECUTABLE)
  find_package(Doxygen)
endif()

if (SPHINX_EXECUTABLE)
  include(blt/cmake/thirdparty/FindSphinx.cmake)
endif()

################################
# Valgrind
################################
include(blt/cmake/thirdparty/FindValgrind.cmake)

################################
# linting via Uncrustify
################################
if (UNCRUSTIFY_EXECUTABLE)
  include(blt/cmake/thirdparty/FindUncrustify.cmake)
endif()

