############################
# Setup compiler options
############################


#############################################
# Support extra compiler flags and defines
#############################################
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)

#
# We don't try to use this approach for CMake generators that support
# multiple configurations. See: CZ JIRA: ATK-45
#
if(NOT CMAKE_CONFIGURATION_TYPES)

    ######################################################
    # Add define we can use when debug builds are enabled
    ######################################################
    if( (CMAKE_BUILD_TYPE MATCHES Debug)
        OR (CMAKE_BUILD_TYPE MATCHES RelWithDebInfo )
      )
        add_definitions(-DATK_DEBUG)
    endif()

    ##########################################
    # Support Extra Flags for the C compiler.
    ##########################################
    if(EXTRA_C_FLAGS)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${EXTRA_C_FLAGS}")
    endif()

    # Extra Flags for the debug builds with the C compiler.
    if(EXTRA_C_FLAGS_DEBUG AND
      ( (CMAKE_BUILD_TYPE MATCHES Debug)
        OR (CMAKE_BUILD_TYPE MATCHES RelWithDebInfo) )
      )
        add_compile_options("${EXTRA_C_FLAGS_DEBUG}")
    endif()

    # Extra Flags for the release builds with the C compiler.
    if(EXTRA_C_FLAGS_RELEASE AND CMAKE_BUILD_TYPE MATCHES RELEASE)
        add_compile_options("${EXTRA_C_FLAGS_RELEASE}")
    endif()

    #############################################
    # Support Extra Flags for the C++ compiler.
    #############################################
    if(EXTRA_CXX_FLAGS)
        add_compile_options("${EXTRA_CXX_FLAGS}")
    endif()

    # Extra Flags for the debug builds with the C++ compiler.
    if(EXTRA_CXX_FLAGS_DEBUG AND
      ( (CMAKE_BUILD_TYPE MATCHES Debug)
        OR (CMAKE_BUILD_TYPE MATCHES RelWithDebInfo) )
      )
        add_compile_options("${EXTRA_CXX_FLAGS_DEBUG}")
    endif()

    # Extra Flags for the release builds with the C++ compiler.
    if(EXTRA_CXX_FLAGS_RELEASE AND CMAKE_BUILD_TYPE MATCHES RELEASE)
        add_compile_options("${EXTRA_CXX_FLAGS_RELEASE}")
    endif()

endif()

################################
# RPath Settings
################################

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
   set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

################################
# Enable C++11 
################################
if (ENABLE_CXX11)
   add_definitions("-DUSE_CXX11")
   
   add_compiler_flag(FLAG_LIST CMAKE_CXX_FLAGS DEFAULT -std=c++11)
   #add_cxx_compiler_flag(-std=c++11)
   #if (NOT HAVE_CXX_FLAG_STD_CXX11)
   #   MESSAGE(FATAL_ERROR "This compiler does not support the standard C++11 flag of '-std=c++11'.  Support for it must be added manually in SetupCompilerOptions.cmake.")
   #endif()
   set(HAVE_CXX_FLAG_STD_CXX11 TRUE)
   MESSAGE(STATUS "Enabled C++11")  
endif()

################################
# Compiler warnings and errors
################################
set(langFlags "CMAKE_C_FLAGS" "CMAKE_CXX_FLAGS")

if (ENABLE_GLOBALCOMPILERWARNINGS)
   MESSAGE(STATUS  "Enabling extra compiler warnings on all targets.")

   foreach(flagVar ${langFlags})   
     add_compiler_flag(FLAG_LIST ${flagVar} 
                     DEFAULT "-W -Wall -Wextra"
                     MSVC "/W4 /Wall /wd4619 /wd4668 /wd4820 /wd4571 /wd4710"
                     )
   endforeach()
endif()

if (ENABLE_GLOBALCOMPILERWARNINGSASERRORS)
   MESSAGE(STATUS  "Enabling treatment of warnings as errors on all targets.")

   foreach(flagVar ${langFlags})   
     add_compiler_flag(FLAG_LIST ${flagVar} 
                     DEFAULT "-Werror"
                     MSVC "/WX"
                     )
   endforeach()
endif()


foreach(flagVar ${langFlags})   
    message(STATUS "${flagVar} flags are:  ${${flagVar}}")
endforeach()


################################
# Enable Fortran
################################
if(ENABLE_FORTRAN)
    add_definitions(-DATK_ENABLE_FORTRAN)

    # if enabled but no fortran compiler, halt the configure
    if(CMAKE_Fortran_COMPILER)
        MESSAGE(STATUS  "Fortran support enabled.")
    else()
        MESSAGE(FATAL_ERROR "Fortran support selected, but no Fortran compiler was found.")
    endif()

    # default property to free form
    set(CMAKE_Fortran_FORMAT FREE)

    # Create macros for Fortran name mangling
    include(FortranCInterface)
    FortranCInterface_HEADER(${HEADER_INCLUDES_DIRECTORY}/common/FC.h MACRO_NAMESPACE "FC_")
else()
    MESSAGE(STATUS  "Fortran support disabled.")
endif()
 

#########################################################################
# Some additional compiler options that can be useful in various targets 
#########################################################################

# Flag for disabling warnings about omp pragmas in the code
add_compiler_flag(FLAG_LIST DISABLE_OMP_PRAGMA_WARNINGS
                  DEFAULT "-Wno-unknown-pragmas"
                  XLC     "-qignprag=omp"
                  )

# Flag for disabling warnings about unused parameters.
# Useful when we include external code.
add_compiler_flag(FLAG_LIST DISABLE_UNUSED_PARAMETER_WARNINGS
                  DEFAULT "-Wno-unused-parameter"
                  XLC     "-qnoinfo=par"
                  )

# message(STATUS "value of DISABLE_OMP_PRAGMA_WARNINGS is ${DISABLE_OMP_PRAGMA_WARNINGS} ")
# message(STATUS "value of DISABLE_UNUSED_PARAMETER_WARNINGS is ${DISABLE_UNUSED_PARAMETER_WARNINGS} ")
 
 
