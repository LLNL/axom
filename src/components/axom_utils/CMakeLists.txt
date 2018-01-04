#------------------------------------------------------------------------------
# Copyright (c) 2017, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-741217
#
# All rights reserved.
#
# This file is part of Axom.
#
# For details about use and distribution, please read axom/LICENSE.
#------------------------------------------------------------------------------

################################
# Add the axom_utils sources
################################
add_subdirectory(src)

################################
# Add tests
################################
if (ENABLE_TESTS)
  add_subdirectory(tests)
endif()

################################
# Add docs
################################
if (ENABLE_DOCS)
  add_subdirectory(docs)
endif()

################################
# Add code checks
################################
axom_add_code_checks(PREFIX axom_utils)

################################
# Create CMake importable
# exports for all of our targets
################################
install(EXPORT axom_utils-targets DESTINATION lib/cmake)

################################
# Compiler checks
################################

if(ENABLE_FORTRAN)

    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/c_loc_with_assumed_shape.f "
! This is expected to fail with gcc 4.7.1, 4.8
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
        ${CMAKE_CURRENT_BINARY_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}/c_loc_with_assumed_shape.f
    )

endif(ENABLE_FORTRAN)
