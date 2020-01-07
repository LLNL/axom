# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# AxomConfig.cmake - Create header of configuration options
#------------------------------------------------------------------------------

## Get Axom version information
message(STATUS "Configuring Axom version ${AXOM_VERSION_FULL}")


## Add a definition to the generated config file for each library dependency
## (optional and built-in) that we might need to know about in the code. We
## check for vars of the form <DEP>_FOUND or ENABLE_<DEP>
set(TPL_DEPS CLI11 CONDUIT CUDA FMT HDF5 MFEM MPI OPENMP RAJA SCR SPARSEHASH UMPIRE )
foreach(dep ${TPL_DEPS})
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


## Add a configuration define for each enabled axom component
foreach(comp ${AXOM_COMPONENTS_ENABLED})
    string(TOUPPER ${comp} comp_uc)
    set(AXOM_USE_${comp_uc} TRUE)
endforeach()

convert_to_native_escaped_file_path(${PROJECT_SOURCE_DIR} AXOM_SRC_DIR)
convert_to_native_escaped_file_path(${CMAKE_BINARY_DIR} AXOM_BIN_DIR)

#------------------------------------------------------------------------------
# Compiler checks
#------------------------------------------------------------------------------

if(ENABLE_FORTRAN)

    file(WRITE ${PROJECT_BINARY_DIR}/c_loc_with_assumed_shape.f "
! This is expected to fail with gcc 4.7.1
! Error: Assumed-shape array 'arg' at (1) cannot be an argument to the
!      procedure 'c_loc' because it is not C interoperable
! https://gcc.gnu.org/bugzilla/show_bug.cgi?id=53945
      program main
      end program main
      subroutine test(arg, addr)
        use iso_c_binding
        integer(C_INT), target, intent(IN) :: arg(:)
        type(C_ptr), intent(OUT) :: addr
        addr = C_LOC(arg)
      end subroutine test
    ")
    try_compile(
        USE_C_LOC_WITH_ASSUMED_SHAPE
        ${PROJECT_BINARY_DIR}
        ${PROJECT_BINARY_DIR}/c_loc_with_assumed_shape.f
    )

endif(ENABLE_FORTRAN)

configure_file(
    ${PROJECT_SOURCE_DIR}/axom/config.hpp.in
    ${CMAKE_BINARY_DIR}/include/axom/config.hpp
)

install(FILES ${CMAKE_BINARY_DIR}/include/axom/config.hpp DESTINATION include/axom)
