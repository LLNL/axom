#Build system author: Chris White (white238@llnl.gov)

#------------------------------------------------------------------------------
# Project setup
#------------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.6)
project(TEMPLATE)
set(PACKAGE_BUGREPORT "hornung@llnl.gov")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules)

#------------------------------------------------------------------------------
# Options
#------------------------------------------------------------------------------
option(BUILD_STATIC_LIBS "Build static libraries" ON)
option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(BUILD_STATIC_EXEC "Build executable linked against static libraries" ON)
option(ENABLE_FORTRAN90 "Build using Fortran 90 components" ON)
option(ENABLE_SUBPROJECT_FOO "Build optional subproject Foo" ON)

#enable C++11
set(CXX_FLAGS ${CXX_FLAGS} "-std=c++11")




#Force atleast static if user turns off both
if(NOT BUILD_STATIC_LIBS AND NOT BUILD_SHARED_LIBS)
  message("Both static and shared libaries were disabled."
    "Building static libraries re-enabled.")
  set(BUILD_STATIC_LIBS ON CACHE BOOL "Build static libraries" FORCE)
endif(NOT BUILD_STATIC_LIBS AND NOT BUILD_SHARED_LIBS)

if(DEFINED ENV{SYS_TYPE})
  set(SYS_TYPE ENV{SYS_TYPE})
else(DEFINED ENV{SYS_TYPE})
  set(SYS_TYPE ${CMAKE_SYSTEM_NAME})
endif(DEFINED ENV{SYS_TYPE})


if(ENABLE_FORTRAN90)
  enable_language(Fortran)
  add_definitions(-DHAVE_F90)
  if("${CMAKE_Fortran_COMPILER}" MATCHES "xlf90")
    add_definitions(-DFTN_ZERO_UNDERSCORE)
  else()
    add_definitions(-DFTN_ONE_UNDERSCORE)
  endif()
endif(ENABLE_FORTRAN90)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#------------------------------------------------------------------------------
# Subprojects
#------------------------------------------------------------------------------
#if(ENABLE_SUBPROJECT_FOO)
#  add_definitions(-DHAVE_FOO)
#  add_subdirectory(foo)
#  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/foo)
#endif(ENABLE_SUBPROJECT_FOO)

#------------------------------------------------------------------------------
# Sources
#------------------------------------------------------------------------------
set(TEMPLATE_HEADERS template.h)

set(TEMPLATE_SOURCES
        template.cpp)

set(MAIN_SOURCES main.cpp)

if(ENABLE_FORTRAN90)
    # add fortrans sources here

    # deal with bgq quirks
    #if("${CMAKE_Fortran_COMPILER}" MATCHES "xlf90")
    #  set_source_files_properties( bar.f PROPERTIES COMPILE_FLAGS -qfixed=132)
    #endif()
  
endif(ENABLE_FORTRAN90)

set(TEMPLATE_INSTALLED_HEADERS
    template.h
    )


#------------------------------------------------------------------------------
# Find Libraries
#------------------------------------------------------------------------------

find_package(hdf5)
if(LIBHDF5_FOUND)
  include_directories(${LIBHDF5_INCLUDE_DIRS})
endif(LIBHDF5_FOUND)

find_package(Lapack)

find_package(Leos)
if(LIBLEOS_FOUND)
  include_directories(${LIBLEOS_INCLUDE_DIRS})
endif(LIBLEOS_FOUND)

#------------------------------------------------------------------------------
# Build Targets
#------------------------------------------------------------------------------
if(BUILD_SHARED_LIBS)
  add_library(template_shared SHARED ${TEMPLATE_SOURCES})
  set_target_properties(template_shared PROPERTIES OUTPUT_NAME template)


  if(LIBHDF5_FOUND)
    target_link_libraries(template_shared ${LIBHDF5_LIBRARIES})
  endif(LIBHDF5_FOUND)

  if(LIBLAPACK_FOUND)
    target_link_libraries(template_shared ${LIBLAPACK_LIBRARIES})
  endif(LIBLAPACK_FOUND)

  if(LIBLEOS_FOUND)
    target_link_libraries(template_shared ${LIBLEOS_LIBRARIES})
  endif(LIBLEOS_FOUND)

endif(BUILD_SHARED_LIBS)

if(BUILD_STATIC_LIBS)
  add_library(template_static STATIC ${TEMPLATE_SOURCES})
  set_target_properties(template_static PROPERTIES OUTPUT_NAME template)

  if(LIBHDF5_FOUND)
    target_link_libraries(template_static ${LIBHDF5_LIBRARIES})
  endif(LIBHDF5_FOUND)

  if(LIBLAPACK_FOUND)
    target_link_libraries(template_static ${LIBLAPACK_LIBRARIES})
  endif(LIBLAPACK_FOUND)

  if(LIBLEOS_FOUND)
    target_link_libraries(template_static ${LIBLEOS_LIBRARIES})
  endif(LIBLEOS_FOUND)

endif(BUILD_STATIC_LIBS)

if (BUILD_STATIC_EXEC) 
    add_executable(main ${MAIN_SOURCES})
    target_link_libraries(main template_static)
endif(BUILD_STATIC_EXEC)
#------------------------------------------------------------------------------
# Testing
#------------------------------------------------------------------------------
if(BUILD_STATIC_LIBS)
 enable_testing()
 add_executable(templatetest ${TEST_SOURCES})
#  add_dependencies(templatetest hdf5)
set_target_properties(templatetest PROPERTIES COMPILE_FLAGS -std=c++11 -lstdc++)
#  target_link_libraries(templatetest template_static)
endif(BUILD_STATIC_LIBS)

#------------------------------------------------------------------------------
# Install files
#------------------------------------------------------------------------------
install(FILES ${TEMPLATE_INSTALLED_HEADERS} DESTINATION include)

if(BUILD_SHARED_LIBS)
  install(TARGETS template_shared DESTINATION lib)
endif(BUILD_SHARED_LIBS)

if(BUILD_STATIC_LIBS)
  install(TARGETS template_static DESTINATION lib)
endif(BUILD_STATIC_LIBS)
