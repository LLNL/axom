! Fortran code that will be inserted into code via shroud splicer blocks

! splicer begin module_use
! map conduit type names to sidre type names
use conduit, only : &
    SIDRE_EMPTY_ID      => CONDUIT_EMPTY_T, &
    SIDRE_INT8_ID       => CONDUIT_INT8_T, &
    SIDRE_INT16_ID      => CONDUIT_INT16_T, &
    SIDRE_INT32_ID      => CONDUIT_INT32_T, &
    SIDRE_INT64_ID      => CONDUIT_INT64_T, &
    SIDRE_UINT8_ID      => CONDUIT_UINT8_T, &
    SIDRE_UINT16_ID     => CONDUIT_UINT16_T, &
    SIDRE_UINT32_ID     => CONDUIT_UINT32_T, &
    SIDRE_UINT64_ID     => CONDUIT_UINT64_T, &
    SIDRE_FLOAT32_ID    => CONDUIT_FLOAT32_T, &
    SIDRE_FLOAT64_ID    => CONDUIT_FLOAT64_T, &
    SIDRE_CHAR8_STR_ID  => CONDUIT_CHAR8_STR_T, &
    SIDRE_INT_ID        => CONDUIT_INT_T, &
    SIDRE_UINT_ID       => CONDUIT_UINT_T, &
    SIDRE_LONG_ID       => CONDUIT_LONG_T, &
    SIDRE_ULONG_ID      => CONDUIT_ULONG_T, &
    SIDRE_FLOAT_ID      => CONDUIT_FLOAT_T, &
    SIDRE_DOUBLE_ID     => CONDUIT_DOUBLE_T
! splicer end module_use

! splicer begin module_top
integer, parameter :: invalid_index = -1
! splicer end module_top


# SIDRE_create_fortran_allocatable_view is not in api.yaml since it is not in src/core and
# only required for the fortran interface.

! splicer begin additional_interfaces
function SIDRE_create_fortran_allocatable_view(group, name, lname, addr, itype, rank) &
   bind(C,name="SIDRE_create_fortran_allocatable_view") &
   result(rv)
      use iso_c_binding
      type(C_PTR), value, intent(IN)    :: group
      character(kind=C_CHAR), intent(IN) :: name(*)
      integer(C_INT), value, intent(IN) :: lname
      type(C_PTR), value                :: addr
      integer(C_INT), value, intent(IN) :: itype
      integer(C_INT), value, intent(IN) :: rank
      type(C_PTR) rv
end function SIDRE_create_fortran_allocatable_view
! splicer end additional_interfaces
