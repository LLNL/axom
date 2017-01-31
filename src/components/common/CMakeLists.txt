################################
# Datastore
################################
project(common)


################################
# Add the common sources
################################
add_subdirectory(src)

################################
# Add tests
################################
if (ENABLE_TESTS)
  add_subdirectory(tests)
endif()


add_code_check_targets(uncrustify.cfg)

################################
# Add docs
################################
if (ENABLE_DOCS)
  add_subdirectory(docs)
endif()

################################
# Create CMake importable
# exports for all of our targets
################################
install(EXPORT ${PROJECT_NAME}-targets DESTINATION lib/cmake)

################################
# Compiler checks
################################

if(ENABLE_FORTRAN)

    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/c_loc_with_assumed_shape.f "
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
        ${CMAKE_CURRENT_BINARY_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}/c_loc_with_assumed_shape.f
    )

endif(ENABLE_FORTRAN)
