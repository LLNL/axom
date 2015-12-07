!
! SidreAllocatablef.f - Routines used by Fortran interface
! Uses cog to insert some generated code into this file.
!

!----------------------------------------------------------------------
![[[cog
!import cog
!import genfsidresplicer as gen
!gen.print_lines(cog.outl, gen.print_sidre_size_allocatable)
!]]]

function sidre_size_allocatable_int_scalar(array) result(rv)
    use iso_c_binding
    implicit none
    integer(C_INT), allocatable, intent(IN) :: array
    integer(C_SIZE_T) :: rv
    rv = 1
end function sidre_size_allocatable_int_scalar

function sidre_size_allocatable_int_1d(array) result(rv)
    use iso_c_binding
    implicit none
    integer(C_INT), allocatable, intent(IN) :: array(:)
    integer(C_SIZE_T) :: rv
    if (allocated(array)) then
        rv = size(array)
    else
        rv = 0
    endif
end function sidre_size_allocatable_int_1d

function sidre_size_allocatable_long_scalar(array) result(rv)
    use iso_c_binding
    implicit none
    integer(C_LONG), allocatable, intent(IN) :: array
    integer(C_SIZE_T) :: rv
    rv = 1
end function sidre_size_allocatable_long_scalar

function sidre_size_allocatable_long_1d(array) result(rv)
    use iso_c_binding
    implicit none
    integer(C_LONG), allocatable, intent(IN) :: array(:)
    integer(C_SIZE_T) :: rv
    if (allocated(array)) then
        rv = size(array)
    else
        rv = 0
    endif
end function sidre_size_allocatable_long_1d

function sidre_size_allocatable_float_scalar(array) result(rv)
    use iso_c_binding
    implicit none
    real(C_FLOAT), allocatable, intent(IN) :: array
    integer(C_SIZE_T) :: rv
    rv = 1
end function sidre_size_allocatable_float_scalar

function sidre_size_allocatable_float_1d(array) result(rv)
    use iso_c_binding
    implicit none
    real(C_FLOAT), allocatable, intent(IN) :: array(:)
    integer(C_SIZE_T) :: rv
    if (allocated(array)) then
        rv = size(array)
    else
        rv = 0
    endif
end function sidre_size_allocatable_float_1d

function sidre_size_allocatable_double_scalar(array) result(rv)
    use iso_c_binding
    implicit none
    real(C_DOUBLE), allocatable, intent(IN) :: array
    integer(C_SIZE_T) :: rv
    rv = 1
end function sidre_size_allocatable_double_scalar

function sidre_size_allocatable_double_1d(array) result(rv)
    use iso_c_binding
    implicit none
    real(C_DOUBLE), allocatable, intent(IN) :: array(:)
    integer(C_SIZE_T) :: rv
    if (allocated(array)) then
        rv = size(array)
    else
        rv = 0
    endif
end function sidre_size_allocatable_double_1d
![[[end]]]

!----------------------------------------------------------------------
![[[cog
!gen.print_lines(cog.outl, gen.print_sidre_address_allocatable)
!]]]

subroutine sidre_address_allocatable_int_scalar(array, addr)
    use iso_c_binding
    implicit none
    integer(C_INT), allocatable, intent(IN), target :: array
    type(C_PTR), intent(OUT) :: addr
    addr = c_loc(array)
end subroutine sidre_address_allocatable_int_scalar

subroutine sidre_address_allocatable_int_1d(array, addr)
    use iso_c_binding
    implicit none
    integer(C_INT), allocatable, intent(IN), target :: array(:)
    type(C_PTR), intent(OUT) :: addr
    addr = c_loc(array)
end subroutine sidre_address_allocatable_int_1d

subroutine sidre_address_allocatable_long_scalar(array, addr)
    use iso_c_binding
    implicit none
    integer(C_LONG), allocatable, intent(IN), target :: array
    type(C_PTR), intent(OUT) :: addr
    addr = c_loc(array)
end subroutine sidre_address_allocatable_long_scalar

subroutine sidre_address_allocatable_long_1d(array, addr)
    use iso_c_binding
    implicit none
    integer(C_LONG), allocatable, intent(IN), target :: array(:)
    type(C_PTR), intent(OUT) :: addr
    addr = c_loc(array)
end subroutine sidre_address_allocatable_long_1d

subroutine sidre_address_allocatable_float_scalar(array, addr)
    use iso_c_binding
    implicit none
    real(C_FLOAT), allocatable, intent(IN), target :: array
    type(C_PTR), intent(OUT) :: addr
    addr = c_loc(array)
end subroutine sidre_address_allocatable_float_scalar

subroutine sidre_address_allocatable_float_1d(array, addr)
    use iso_c_binding
    implicit none
    real(C_FLOAT), allocatable, intent(IN), target :: array(:)
    type(C_PTR), intent(OUT) :: addr
    addr = c_loc(array)
end subroutine sidre_address_allocatable_float_1d

subroutine sidre_address_allocatable_double_scalar(array, addr)
    use iso_c_binding
    implicit none
    real(C_DOUBLE), allocatable, intent(IN), target :: array
    type(C_PTR), intent(OUT) :: addr
    addr = c_loc(array)
end subroutine sidre_address_allocatable_double_scalar

subroutine sidre_address_allocatable_double_1d(array, addr)
    use iso_c_binding
    implicit none
    real(C_DOUBLE), allocatable, intent(IN), target :: array(:)
    type(C_PTR), intent(OUT) :: addr
    addr = c_loc(array)
end subroutine sidre_address_allocatable_double_1d
![[[end]]]

!----------------------------------------------------------------------
![[[cog
!gen.print_lines(cog.outl, gen.print_sidre_allocate_allocatable)
!]]]

subroutine sidre_allocate_allocatable_int_scalar(array, nitems)
    use iso_c_binding
    implicit none
    integer(C_INT), allocatable, intent(OUT), target :: array
    integer(C_INT), intent(IN) :: nitems
    allocate(array)
end subroutine sidre_allocate_allocatable_int_scalar

subroutine sidre_allocate_allocatable_int_1d(array, nitems)
    use iso_c_binding
    implicit none
    integer(C_INT), allocatable, intent(OUT), target :: array(:)
    integer(C_INT), intent(IN) :: nitems
    allocate(array(nitems))
end subroutine sidre_allocate_allocatable_int_1d

subroutine sidre_allocate_allocatable_long_scalar(array, nitems)
    use iso_c_binding
    implicit none
    integer(C_LONG), allocatable, intent(OUT), target :: array
    integer(C_INT), intent(IN) :: nitems
    allocate(array)
end subroutine sidre_allocate_allocatable_long_scalar

subroutine sidre_allocate_allocatable_long_1d(array, nitems)
    use iso_c_binding
    implicit none
    integer(C_LONG), allocatable, intent(OUT), target :: array(:)
    integer(C_INT), intent(IN) :: nitems
    allocate(array(nitems))
end subroutine sidre_allocate_allocatable_long_1d

subroutine sidre_allocate_allocatable_float_scalar(array, nitems)
    use iso_c_binding
    implicit none
    real(C_FLOAT), allocatable, intent(OUT), target :: array
    integer(C_INT), intent(IN) :: nitems
    allocate(array)
end subroutine sidre_allocate_allocatable_float_scalar

subroutine sidre_allocate_allocatable_float_1d(array, nitems)
    use iso_c_binding
    implicit none
    real(C_FLOAT), allocatable, intent(OUT), target :: array(:)
    integer(C_INT), intent(IN) :: nitems
    allocate(array(nitems))
end subroutine sidre_allocate_allocatable_float_1d

subroutine sidre_allocate_allocatable_double_scalar(array, nitems)
    use iso_c_binding
    implicit none
    real(C_DOUBLE), allocatable, intent(OUT), target :: array
    integer(C_INT), intent(IN) :: nitems
    allocate(array)
end subroutine sidre_allocate_allocatable_double_scalar

subroutine sidre_allocate_allocatable_double_1d(array, nitems)
    use iso_c_binding
    implicit none
    real(C_DOUBLE), allocatable, intent(OUT), target :: array(:)
    integer(C_INT), intent(IN) :: nitems
    allocate(array(nitems))
end subroutine sidre_allocate_allocatable_double_1d
![[[end]]]

!----------------------------------------------------------------------
![[[cog
!gen.print_lines(cog.outl, gen.print_sidre_deallocate_allocatable)
!]]]

subroutine sidre_deallocate_allocatable_int_scalar(array)
    use iso_c_binding
    implicit none
    integer(C_INT), allocatable, intent(INOUT), target :: array
    deallocate(array)
end subroutine sidre_deallocate_allocatable_int_scalar

subroutine sidre_deallocate_allocatable_int_1d(array)
    use iso_c_binding
    implicit none
    integer(C_INT), allocatable, intent(INOUT), target :: array(:)
    deallocate(array)
end subroutine sidre_deallocate_allocatable_int_1d

subroutine sidre_deallocate_allocatable_long_scalar(array)
    use iso_c_binding
    implicit none
    integer(C_LONG), allocatable, intent(INOUT), target :: array
    deallocate(array)
end subroutine sidre_deallocate_allocatable_long_scalar

subroutine sidre_deallocate_allocatable_long_1d(array)
    use iso_c_binding
    implicit none
    integer(C_LONG), allocatable, intent(INOUT), target :: array(:)
    deallocate(array)
end subroutine sidre_deallocate_allocatable_long_1d

subroutine sidre_deallocate_allocatable_float_scalar(array)
    use iso_c_binding
    implicit none
    real(C_FLOAT), allocatable, intent(INOUT), target :: array
    deallocate(array)
end subroutine sidre_deallocate_allocatable_float_scalar

subroutine sidre_deallocate_allocatable_float_1d(array)
    use iso_c_binding
    implicit none
    real(C_FLOAT), allocatable, intent(INOUT), target :: array(:)
    deallocate(array)
end subroutine sidre_deallocate_allocatable_float_1d

subroutine sidre_deallocate_allocatable_double_scalar(array)
    use iso_c_binding
    implicit none
    real(C_DOUBLE), allocatable, intent(INOUT), target :: array
    deallocate(array)
end subroutine sidre_deallocate_allocatable_double_scalar

subroutine sidre_deallocate_allocatable_double_1d(array)
    use iso_c_binding
    implicit none
    real(C_DOUBLE), allocatable, intent(INOUT), target :: array(:)
    deallocate(array)
end subroutine sidre_deallocate_allocatable_double_1d
![[[end]]]

!----------------------------------------------------------------------
![[[cog
!gen.print_lines(cog.outl, gen.print_sidre_reallocate_allocatable)
!]]]

subroutine sidre_reallocate_allocatable_int_scalar(array, nitems)
    use iso_c_binding
    implicit none
    integer(C_INT), allocatable, intent(INOUT), target :: array
    integer(C_INT), value, intent(IN) :: nitems
    allocate(array)
! XXX move_alloc ...
end subroutine sidre_reallocate_allocatable_int_scalar

subroutine sidre_reallocate_allocatable_int_1d(array, nitems)
    use iso_c_binding
    implicit none
    integer(C_INT), allocatable, intent(INOUT), target :: array(:)
    integer(C_INT), value, intent(IN) :: nitems
    allocate(array(nitems))
! XXX move_alloc ...
end subroutine sidre_reallocate_allocatable_int_1d

subroutine sidre_reallocate_allocatable_long_scalar(array, nitems)
    use iso_c_binding
    implicit none
    integer(C_LONG), allocatable, intent(INOUT), target :: array
    integer(C_INT), value, intent(IN) :: nitems
    allocate(array)
! XXX move_alloc ...
end subroutine sidre_reallocate_allocatable_long_scalar

subroutine sidre_reallocate_allocatable_long_1d(array, nitems)
    use iso_c_binding
    implicit none
    integer(C_LONG), allocatable, intent(INOUT), target :: array(:)
    integer(C_INT), value, intent(IN) :: nitems
    allocate(array(nitems))
! XXX move_alloc ...
end subroutine sidre_reallocate_allocatable_long_1d

subroutine sidre_reallocate_allocatable_float_scalar(array, nitems)
    use iso_c_binding
    implicit none
    real(C_FLOAT), allocatable, intent(INOUT), target :: array
    integer(C_INT), value, intent(IN) :: nitems
    allocate(array)
! XXX move_alloc ...
end subroutine sidre_reallocate_allocatable_float_scalar

subroutine sidre_reallocate_allocatable_float_1d(array, nitems)
    use iso_c_binding
    implicit none
    real(C_FLOAT), allocatable, intent(INOUT), target :: array(:)
    integer(C_INT), value, intent(IN) :: nitems
    allocate(array(nitems))
! XXX move_alloc ...
end subroutine sidre_reallocate_allocatable_float_1d

subroutine sidre_reallocate_allocatable_double_scalar(array, nitems)
    use iso_c_binding
    implicit none
    real(C_DOUBLE), allocatable, intent(INOUT), target :: array
    integer(C_INT), value, intent(IN) :: nitems
    allocate(array)
! XXX move_alloc ...
end subroutine sidre_reallocate_allocatable_double_scalar

subroutine sidre_reallocate_allocatable_double_1d(array, nitems)
    use iso_c_binding
    implicit none
    real(C_DOUBLE), allocatable, intent(INOUT), target :: array(:)
    integer(C_INT), value, intent(IN) :: nitems
    allocate(array(nitems))
! XXX move_alloc ...
end subroutine sidre_reallocate_allocatable_double_1d
![[[end]]]

