####################################
# Datastore 3rd Party Dependencies
####################################

################################
# Conduit
################################
include(cmake/FindConduit.cmake)


################################
# Sparsehash
################################
include(cmake/FindSparsehash.cmake)

################################
# Documentation Packages
################################

find_package(Doxygen)
include(cmake/FindSphinx.cmake)

################################
# linting via Uncrustify
################################
include(cmake/FindUncrustify.cmake)

################################
# Find boost headers
################################
if (ENABLE_BOOST)
  if (DEFINED BOOST_ROOT)
    find_package(Boost
                 1.55
                 REQUIRED)
    MESSAGE(STATUS "Boost include dir: " ${Boost_INCLUDE_DIR})
    MESSAGE(STATUS "Boost version: " ${Boost_VERSION} )
  else()
    MESSAGE(FATAL_ERROR "ENABLE_BOOST is true, but BOOST_ROOT was not set.  Check your host-config file.")
  endif()
endif()

################################
# Python
################################

if(PYTHON_EXECUTABLE)
    ################################
    # Setup includes for Python & Numpy
    ################################
    include(FindPython)
    message(STATUS "Using Python Include: ${PYTHON_INCLUDE_DIRS}")
    include_directories(${PYTHON_INCLUDE_DIRS})
    # if we don't find python, throw a fatal error
    if(NOT PYTHON_FOUND)
        message(FATAL_ERROR "ENABLE_EXECUTABLE is set, but Python wasn't found.")
    endif()

    ## Set the Python module directory
    # relative path (used with install)
    set(CMAKE_Python_SITE_PACKAGES
        "lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages"
        CACHE PATH
        "Relative path where all Python modules will go in the install tree"
    )
    # build site-packages
    set(CMAKE_Python_MODULE_DIRECTORY
        "${PROJECT_BINARY_DIR}/${CMAKE_Python_SITE_PACKAGES}"
        CACHE PATH
        "Directory where all Python modules will go in the build tree"
    )

    file(MAKE_DIRECTORY ${CMAKE_Python_MODULE_DIRECTORY})
    set(ENV{PYTHONPATH} ${CMAKE_Python_MODULE_DIRECTORY})

    INSTALL(DIRECTORY DESTINATION ${CMAKE_Python_SITE_PACKAGES})
    INSTALL(CODE " set(ENV\{PYTHONPATH\} ${CMAKE_INSTALL_PREFIX}/${CMAKE_Python_SITE_PACKAGES}) ")
endif(PYTHON_EXECUTABLE)
