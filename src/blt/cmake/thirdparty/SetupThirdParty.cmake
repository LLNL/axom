####################################
# BLT 3rd Party Lib Support
####################################

################################
# MPI
################################
message(STATUS "MPI Support is ${ENABLE_MPI}")
if (ENABLE_MPI)
  find_package(MPI REQUIRED)
  message(STATUS "MPI C Compile Flags: ${MPI_C_COMPILE_FLAGS}")
  message(STATUS "MPI C Include Path:  ${MPI_C_INCLUDE_PATH}")
  message(STATUS "MPI C Link Flags:    ${MPI_C_LINK_FLAGS}")
  message(STATUS "MPI C Libraries:     ${MPI_C_LIBRARIES}")
endif()

################################
# OpenMP
################################
message(STATUS "OpenMP Support is ${ENABLE_OPENMP}")
if(ENABLE_OPENMP)
    find_package(OpenMP REQUIRED)
    message(STATUS "OpenMP CXX Flags: ${OpenMP_CXX_FLAGS}")
endif()


################################
# Documentation Packages
################################
if (DOXYGEN_EXECUTABLE)
  find_package(Doxygen)
endif()

if (SPHINX_EXECUTABLE)
  include(blt/cmake/thirdparty/FindSphinx.cmake)
endif()


if (SPHINX_EXECUTABLE)
  include(blt/cmake/thirdparty/FindValgrind.cmake)
endif()


################################
# linting via Uncrustify
################################
if (UNCRUSTIFY_EXECUTABLE)
  include(blt/cmake/thirdparty/FindUncrustify.cmake)
endif()

