!
! SidreAllocatablef.f - Routines used by Fortran interface
! Uses cog to insert some generated code into this file.
!

!----------------------------------------------------------------------
![[[cog
!import cog
!import genfsidresplicer as gen
!gen.print_lines(cog.outl, gen.print_atk_size_allocatable)
!]]]

function atk_size_allocatable_int_scalar_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    integer(C_INT), allocatable, intent(IN) :: array
    integer(C_SIZE_T) :: rv
    rv = 1
end function atk_size_allocatable_int_scalar_ptr

function atk_size_allocatable_int_1d_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    integer(C_INT), allocatable, intent(IN) :: array(:)
    integer(C_SIZE_T) :: rv
    if (allocated(array)) then
        rv = size(array)
    else
        rv = 0
    endif
end function atk_size_allocatable_int_1d_ptr

function atk_size_allocatable_long_scalar_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    integer(C_LONG), allocatable, intent(IN) :: array
    integer(C_SIZE_T) :: rv
    rv = 1
end function atk_size_allocatable_long_scalar_ptr

function atk_size_allocatable_long_1d_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    integer(C_LONG), allocatable, intent(IN) :: array(:)
    integer(C_SIZE_T) :: rv
    if (allocated(array)) then
        rv = size(array)
    else
        rv = 0
    endif
end function atk_size_allocatable_long_1d_ptr

function atk_size_allocatable_float_scalar_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    real(C_FLOAT), allocatable, intent(IN) :: array
    integer(C_SIZE_T) :: rv
    rv = 1
end function atk_size_allocatable_float_scalar_ptr

function atk_size_allocatable_float_1d_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    real(C_FLOAT), allocatable, intent(IN) :: array(:)
    integer(C_SIZE_T) :: rv
    if (allocated(array)) then
        rv = size(array)
    else
        rv = 0
    endif
end function atk_size_allocatable_float_1d_ptr

function atk_size_allocatable_double_scalar_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    real(C_DOUBLE), allocatable, intent(IN) :: array
    integer(C_SIZE_T) :: rv
    rv = 1
end function atk_size_allocatable_double_scalar_ptr

function atk_size_allocatable_double_1d_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    real(C_DOUBLE), allocatable, intent(IN) :: array(:)
    integer(C_SIZE_T) :: rv
    if (allocated(array)) then
        rv = size(array)
    else
        rv = 0
    endif
end function atk_size_allocatable_double_1d_ptr
![[[end]]]

!----------------------------------------------------------------------
![[[cog
!gen.print_lines(cog.outl, gen.print_atk_address_allocatable)
!]]]

function atk_address_allocatable_int_scalar_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    integer(C_INT), allocatable, intent(IN), target :: array
    type(C_PTR) :: rv
    rv = c_loc(array)
end function atk_address_allocatable_int_scalar_ptr

function atk_address_allocatable_int_1d_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    integer(C_INT), allocatable, intent(IN), target :: array(:)
    type(C_PTR) :: rv
    rv = c_loc(array)
end function atk_address_allocatable_int_1d_ptr

function atk_address_allocatable_long_scalar_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    integer(C_LONG), allocatable, intent(IN), target :: array
    type(C_PTR) :: rv
    rv = c_loc(array)
end function atk_address_allocatable_long_scalar_ptr

function atk_address_allocatable_long_1d_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    integer(C_LONG), allocatable, intent(IN), target :: array(:)
    type(C_PTR) :: rv
    rv = c_loc(array)
end function atk_address_allocatable_long_1d_ptr

function atk_address_allocatable_float_scalar_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    real(C_FLOAT), allocatable, intent(IN), target :: array
    type(C_PTR) :: rv
    rv = c_loc(array)
end function atk_address_allocatable_float_scalar_ptr

function atk_address_allocatable_float_1d_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    real(C_FLOAT), allocatable, intent(IN), target :: array(:)
    type(C_PTR) :: rv
    rv = c_loc(array)
end function atk_address_allocatable_float_1d_ptr

function atk_address_allocatable_double_scalar_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    real(C_DOUBLE), allocatable, intent(IN), target :: array
    type(C_PTR) :: rv
    rv = c_loc(array)
end function atk_address_allocatable_double_scalar_ptr

function atk_address_allocatable_double_1d_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    real(C_DOUBLE), allocatable, intent(IN), target :: array(:)
    type(C_PTR) :: rv
    rv = c_loc(array)
end function atk_address_allocatable_double_1d_ptr
![[[end]]]

!----------------------------------------------------------------------
![[[cog
!gen.print_lines(cog.outl, gen.print_atk_allocate_allocatable)
!]]]

subroutine atk_allocate_allocatable_int_scalar_ptr(array, nitems)
    use iso_c_binding
    implicit none
    integer(C_INT), allocatable, intent(OUT), target :: array
    integer(C_INT), value, intent(IN) :: nitems
    allocate(array)
end subroutine atk_allocate_allocatable_int_scalar_ptr

subroutine atk_allocate_allocatable_int_1d_ptr(array, nitems)
    use iso_c_binding
    implicit none
    integer(C_INT), allocatable, intent(OUT), target :: array(:)
    integer(C_INT), value, intent(IN) :: nitems
    allocate(array(nitems))
end subroutine atk_allocate_allocatable_int_1d_ptr

subroutine atk_allocate_allocatable_long_scalar_ptr(array, nitems)
    use iso_c_binding
    implicit none
    integer(C_LONG), allocatable, intent(OUT), target :: array
    integer(C_INT), value, intent(IN) :: nitems
    allocate(array)
end subroutine atk_allocate_allocatable_long_scalar_ptr

subroutine atk_allocate_allocatable_long_1d_ptr(array, nitems)
    use iso_c_binding
    implicit none
    integer(C_LONG), allocatable, intent(OUT), target :: array(:)
    integer(C_INT), value, intent(IN) :: nitems
    allocate(array(nitems))
end subroutine atk_allocate_allocatable_long_1d_ptr

subroutine atk_allocate_allocatable_float_scalar_ptr(array, nitems)
    use iso_c_binding
    implicit none
    real(C_FLOAT), allocatable, intent(OUT), target :: array
    integer(C_INT), value, intent(IN) :: nitems
    allocate(array)
end subroutine atk_allocate_allocatable_float_scalar_ptr

subroutine atk_allocate_allocatable_float_1d_ptr(array, nitems)
    use iso_c_binding
    implicit none
    real(C_FLOAT), allocatable, intent(OUT), target :: array(:)
    integer(C_INT), value, intent(IN) :: nitems
    allocate(array(nitems))
end subroutine atk_allocate_allocatable_float_1d_ptr

subroutine atk_allocate_allocatable_double_scalar_ptr(array, nitems)
    use iso_c_binding
    implicit none
    real(C_DOUBLE), allocatable, intent(OUT), target :: array
    integer(C_INT), value, intent(IN) :: nitems
    allocate(array)
end subroutine atk_allocate_allocatable_double_scalar_ptr

subroutine atk_allocate_allocatable_double_1d_ptr(array, nitems)
    use iso_c_binding
    implicit none
    real(C_DOUBLE), allocatable, intent(OUT), target :: array(:)
    integer(C_INT), value, intent(IN) :: nitems
    allocate(array(nitems))
end subroutine atk_allocate_allocatable_double_1d_ptr
![[[end]]]

!----------------------------------------------------------------------
![[[cog
!gen.print_lines(cog.outl, gen.print_atk_deallocate_allocatable)
!]]]

subroutine atk_deallocate_allocatable_int_scalar_ptr(array)
    use iso_c_binding
    implicit none
    integer(C_INT), allocatable, intent(INOUT), target :: array
    deallocate(array)
end subroutine atk_deallocate_allocatable_int_scalar_ptr

subroutine atk_deallocate_allocatable_int_1d_ptr(array)
    use iso_c_binding
    implicit none
    integer(C_INT), allocatable, intent(INOUT), target :: array(:)
    deallocate(array)
end subroutine atk_deallocate_allocatable_int_1d_ptr

subroutine atk_deallocate_allocatable_long_scalar_ptr(array)
    use iso_c_binding
    implicit none
    integer(C_LONG), allocatable, intent(INOUT), target :: array
    deallocate(array)
end subroutine atk_deallocate_allocatable_long_scalar_ptr

subroutine atk_deallocate_allocatable_long_1d_ptr(array)
    use iso_c_binding
    implicit none
    integer(C_LONG), allocatable, intent(INOUT), target :: array(:)
    deallocate(array)
end subroutine atk_deallocate_allocatable_long_1d_ptr

subroutine atk_deallocate_allocatable_float_scalar_ptr(array)
    use iso_c_binding
    implicit none
    real(C_FLOAT), allocatable, intent(INOUT), target :: array
    deallocate(array)
end subroutine atk_deallocate_allocatable_float_scalar_ptr

subroutine atk_deallocate_allocatable_float_1d_ptr(array)
    use iso_c_binding
    implicit none
    real(C_FLOAT), allocatable, intent(INOUT), target :: array(:)
    deallocate(array)
end subroutine atk_deallocate_allocatable_float_1d_ptr

subroutine atk_deallocate_allocatable_double_scalar_ptr(array)
    use iso_c_binding
    implicit none
    real(C_DOUBLE), allocatable, intent(INOUT), target :: array
    deallocate(array)
end subroutine atk_deallocate_allocatable_double_scalar_ptr

subroutine atk_deallocate_allocatable_double_1d_ptr(array)
    use iso_c_binding
    implicit none
    real(C_DOUBLE), allocatable, intent(INOUT), target :: array(:)
    deallocate(array)
end subroutine atk_deallocate_allocatable_double_1d_ptr
![[[end]]]

!----------------------------------------------------------------------
![[[cog
!gen.print_lines(cog.outl, gen.print_atk_reallocate_allocatable)
!]]]

subroutine atk_reallocate_allocatable_int_scalar_ptr(array, nitems)
    use iso_c_binding
    implicit none
    integer(C_INT), allocatable, intent(INOUT), target :: array
    integer(C_INT), value, intent(IN) :: nitems
    allocate(array)
! XXX move_alloc ...
end subroutine atk_reallocate_allocatable_int_scalar_ptr

subroutine atk_reallocate_allocatable_int_1d_ptr(array, nitems)
    use iso_c_binding
    implicit none
    integer(C_INT), allocatable, intent(INOUT), target :: array(:)
    integer(C_INT), value, intent(IN) :: nitems
    allocate(array(nitems))
! XXX move_alloc ...
end subroutine atk_reallocate_allocatable_int_1d_ptr

subroutine atk_reallocate_allocatable_long_scalar_ptr(array, nitems)
    use iso_c_binding
    implicit none
    integer(C_LONG), allocatable, intent(INOUT), target :: array
    integer(C_INT), value, intent(IN) :: nitems
    allocate(array)
! XXX move_alloc ...
end subroutine atk_reallocate_allocatable_long_scalar_ptr

subroutine atk_reallocate_allocatable_long_1d_ptr(array, nitems)
    use iso_c_binding
    implicit none
    integer(C_LONG), allocatable, intent(INOUT), target :: array(:)
    integer(C_INT), value, intent(IN) :: nitems
    allocate(array(nitems))
! XXX move_alloc ...
end subroutine atk_reallocate_allocatable_long_1d_ptr

subroutine atk_reallocate_allocatable_float_scalar_ptr(array, nitems)
    use iso_c_binding
    implicit none
    real(C_FLOAT), allocatable, intent(INOUT), target :: array
    integer(C_INT), value, intent(IN) :: nitems
    allocate(array)
! XXX move_alloc ...
end subroutine atk_reallocate_allocatable_float_scalar_ptr

subroutine atk_reallocate_allocatable_float_1d_ptr(array, nitems)
    use iso_c_binding
    implicit none
    real(C_FLOAT), allocatable, intent(INOUT), target :: array(:)
    integer(C_INT), value, intent(IN) :: nitems
    allocate(array(nitems))
! XXX move_alloc ...
end subroutine atk_reallocate_allocatable_float_1d_ptr

subroutine atk_reallocate_allocatable_double_scalar_ptr(array, nitems)
    use iso_c_binding
    implicit none
    real(C_DOUBLE), allocatable, intent(INOUT), target :: array
    integer(C_INT), value, intent(IN) :: nitems
    allocate(array)
! XXX move_alloc ...
end subroutine atk_reallocate_allocatable_double_scalar_ptr

subroutine atk_reallocate_allocatable_double_1d_ptr(array, nitems)
    use iso_c_binding
    implicit none
    real(C_DOUBLE), allocatable, intent(INOUT), target :: array(:)
    integer(C_INT), value, intent(IN) :: nitems
    allocate(array(nitems))
! XXX move_alloc ...
end subroutine atk_reallocate_allocatable_double_1d_ptr
![[[end]]]

