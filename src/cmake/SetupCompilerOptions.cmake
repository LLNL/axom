############################
# Setup compiler options
############################

#####################################################
# Set some variables to simplify determining compiler
# Compiler string list from: 
#   https://cmake.org/cmake/help/v3.0/variable/CMAKE_LANG_COMPILER_ID.html
####################################################3

if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU") 
    set(COMPILER_FAMILY_IS_GNU 1)
    message(STATUS "Compiler family is GNU")
    
elseif(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang") # For Clang or AppleClang
    set(COMPILER_FAMILY_IS_CLANG 1)
    message(STATUS "Compiler family is Clang")
    
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "XL") 
    set(COMPILER_FAMILY_IS_XL 1)
    message(STATUS "Compiler family is XL")    
    
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel") 
    set(COMPILER_FAMILY_IS_INTEL 1)
    message(STATUS "Compiler family is Intel")
   
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "MSVC") 
    set(COMPILER_FAMILY_IS_MSVC 1)
    message(STATUS "Compiler family is MSVC")
    
endif()




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
   
   blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -std=c++11)
   set(HAVE_CXX_FLAG_STD_CXX11 TRUE)
   
   MESSAGE(STATUS "C++11 support is ON")  
else()
   MESSAGE(STATUS "C++11 support is OFF")  
endif()

##################################################################
# Additional compiler warnings and treatment of warnings as errors
##################################################################


blt_append_custom_compiler_flag(
    FLAGS_VAR ATK_ENABLE_ALL_WARNINGS
     DEFAULT "-Wall -Wextra"
     CLANG   "-Wall -Wextra" 
                    # Additional  possibilities for clang include: 
                    #       "-Wdocumentation -Wdeprecated -Weverything"
     MSVC    "/W4 /Wall /wd4619 /wd4668 /wd4820 /wd4571 /wd4710"
     XL      ""     # qinfo=<grp> produces additional messages on XL
                    # qflag=<x>:<x> defines min severity level to produce messages on XL
                    #     where x is i info, w warning, e error, s severe; default is: 
                    # (default is  qflag=i:i)
     )

blt_append_custom_compiler_flag(
    FLAGS_VAR ATK_TREAT_WARNINGS_AS_ERRORS
     DEFAULT  "-Werror"
     MSVC     "/WX"
     XL       "qhalt=w"       # i info, w warning, e error, s severe (default)
     )

set(langFlags "CMAKE_C_FLAGS" "CMAKE_CXX_FLAGS")

if (ENABLE_ALL_WARNINGS)
   MESSAGE(STATUS  "Enabling all compiler warnings on all targets.")

   foreach(flagVar ${langFlags})
     set(${flagVar} "${${flagVar}} ${ATK_ENABLE_ALL_WARNINGS}") 
   endforeach()
endif()

if (ENABLE_WARNINGS_AS_ERRORS)
   MESSAGE(STATUS  "Enabling treatment of warnings as errors on all targets.")

   foreach(flagVar ${langFlags})   
     set(${flagVar} "${${flagVar}} ${ATK_TREAT_WARNINGS_AS_ERRORS}") 
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

    if (ENABLE_MPI)
        # Determine if we should use fortran mpif.h header or fortran mpi module
        find_path(mpif_path
            NAMES "mpif.h"
            PATHS ${MPI_Fortran_INCLUDE_PATH}
            NO_DEFAULT_PATH
            )
        
        if(mpif_path)
            set(MPI_Fortran_USE_MPIF ON CACHE PATH "")
            message(STATUS "Using MPI Fortran header: mpif.h")
        else()
            set(MPI_Fortran_USE_MPIF OFF CACHE PATH "")
            message(STATUS "Using MPI Fortran module: mpi.mod")
        endif()
    endif()
    
else()
    MESSAGE(STATUS  "Fortran support disabled.")
endif()
 

##############################################################################
# Setup some additional compiler options that can be useful in various targets
# These are stored in their own variables.
# Usage: To add one of these sets of flags to some source files:
#   get_source_file_property(_origflags <src_file> COMPILE_FLAGS)
#   set_source_files_properties(<list_of_src_files> 
#        PROPERTIES COMPILE_FLAGS "${_origFlags} ${<flags_variable}" )
##############################################################################

set(custom_compiler_flags_list) # Tracks custom compiler flags for logging

# Flag for disabling warnings about omp pragmas in the code
blt_append_custom_compiler_flag(FLAGS_VAR ATK_DISABLE_OMP_PRAGMA_WARNINGS
                  DEFAULT "-Wno-unknown-pragmas"
                  XL      "-qignprag=omp"
                  INTEL   "-diag-disable 3180"
                  )
list(APPEND custom_compiler_flags_list ATK_DISABLE_OMP_PRAGMA_WARNINGS)

# Flag for disabling warnings about unused parameters.
# Useful when we include external code.
blt_append_custom_compiler_flag(FLAGS_VAR ATK_DISABLE_UNUSED_PARAMETER_WARNINGS
                  DEFAULT "-Wno-unused-parameter"
                  XL      "-qinfo=nopar"
                  )
list(APPEND custom_compiler_flags_list ATK_DISABLE_UNUSED_PARAMETER_WARNINGS)

# Flag for disabling warnings about unused variables
# Useful when we include external code.
blt_append_custom_compiler_flag(FLAGS_VAR ATK_DISABLE_UNUSED_VARIABLE_WARNINGS
                  DEFAULT "-Wno-unused-variable"
                  XL      "-qinfo=nouse"
                  )
list(APPEND custom_compiler_flags_list ATK_DISABLE_UNUSED_VARIABLE_WARNINGS)

# Flag for disabling warnings about variables that may be uninitialized.
# Useful when we are using compiler generated interface code (e.g. in shroud)
blt_append_custom_compiler_flag(FLAGS_VAR ATK_DISABLE_UNINITIALIZED_WARNINGS
                  DEFAULT "-Wno-uninitialized"
                  XL      "-qsuppress=1540-1102"
                  )
list(APPEND custom_compiler_flags_list ATK_DISABLE_UNINITIALIZED_WARNINGS)

# Flag for disabling warnings about strict aliasing.
# Useful when we are using compiler generated interface code (e.g. in shroud)
blt_append_custom_compiler_flag(FLAGS_VAR ATK_DISABLE_ALIASING_WARNINGS
                  DEFAULT "-Wno-strict-aliasing"
                  XL      ""
                  )
list(APPEND custom_compiler_flags_list ATK_DISABLE_ALIASING_WARNINGS)

# Flag for enabling the C preprocessor in fortran.
# (Note KW 5/2016) The XL flag only applies to *.f90 files -- I could not find a more general solution.   
# 		xlf only allows one file remapping at a time. If you have *.f files, '-qsuffix=cpp=f' should work.
# 		Alternatively, you can rename the file's extension to automatically invoke the preprocessor (e.g. *.f90 ->  *.F90)
blt_append_custom_compiler_flag(FLAGS_VAR ATK_PREPROCESS_FORTRAN
                  DEFAULT "-cpp"
                  XL      "-qsuffix=cpp=f90" 	# Note: only invokes the preprocessor on files with extension *.f90
                  )
list(APPEND custom_compiler_flags_list ATK_PREPROCESS_FORTRAN)

   
# message(STATUS "Custom compiler flags:")
# foreach(flag ${custom_compiler_flags_list})
#    message(STATUS "\tvalue of ${flag} is '${${flag}}'")
# endforeach()
 
