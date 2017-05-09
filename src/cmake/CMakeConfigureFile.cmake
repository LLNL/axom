#
# CMakeConfigureFile.cmake - Create header of configuration options
#

if( (CMAKE_CXX_STANDARD EQUAL 11) OR (CMAKE_CXX_STANDARD EQUAL 14) )
    set(AXOM_USE_CXX11 TRUE)
endif()


## Add a definition to the generated config file for each library dependency
## (optional and built-in) that we might need to know about in the code
## Note: BLT adds USE_MPI and USE_OPENMP as compile define flags for targets
##       that are configured with MPI and OPENMP, respectively.

set(TPL_DEPS CONDUIT HDF5 SPARSEHASH FMT BOOST MPI MFEM)  # vars of the form DEP_FOUND
foreach(dep in ${TPL_DEPS})
    if( ${dep}_FOUND OR ENABLE_${dep} )
        set(AXOM_USE_${dep} TRUE  )
    endif()
endforeach()

# Handle MPI Fortran headers
if(ENABLE_MPI AND ENABLE_FORTRAN)
  if(MPI_Fortran_USE_MPIF)
    set(AXOM_USE_MPIF_HEADER TRUE)
  endif()
endif()


# If Sparsehash was found, AXOM_USE_SPARSEHASH was set above in the TPL_DEPS
# loop.  If not, we must use a standard container--std::unordered_map when
# using C++11, std::map otherwise.  std::map is expected to perform poorly
# with large amounts of Sidre objects, so it is recommended to make sure
# Sparsehash is available for non-C++ 11 builds.
if(NOT AXOM_USE_SPARSEHASH)
  if(AXOM_USE_CXX11)
    set(AXOM_USE_STD_UNORDERED_MAP TRUE)
  endif()
endif()


## Add a configuration define for each enabled axom component
set(COMPS AXOM_UTILS LUMBERJACK SLIC SLAM SIDRE MINT PRIMAL QUEST SPIO)
foreach(comp in ${COMPS})
    if( ENABLE_${comp} )
        set(AXOM_USE_${comp} TRUE)
    endif()
endforeach()



configure_file(
    include/config.hpp.in
    ${HEADER_INCLUDES_DIRECTORY}/axom/config.hpp
)
