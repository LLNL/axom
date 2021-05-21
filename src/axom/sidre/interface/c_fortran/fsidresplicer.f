! Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
! other Axom Project Developers. See the top-level LICENSE file for details.
!
! SPDX-License-Identifier: (BSD-3-Clause)

! Fortran code that will be inserted into code via shroud splicer blocks

! splicer begin file_top
#include "axom/config.hpp"
! splicer end file_top

! splicer begin module_use
! map conduit type names to sidre type names
use conduit, only : &
    SIDRE_NO_TYPE_ID    => CONDUIT_EMPTY_ID, &
    SIDRE_INT8_ID       => CONDUIT_INT8_ID, &
    SIDRE_INT16_ID      => CONDUIT_INT16_ID, &
    SIDRE_INT32_ID      => CONDUIT_INT32_ID, &
    SIDRE_INT64_ID      => CONDUIT_INT64_ID, &
    SIDRE_UINT8_ID      => CONDUIT_UINT8_ID, &
    SIDRE_UINT16_ID     => CONDUIT_UINT16_ID, &
    SIDRE_UINT32_ID     => CONDUIT_UINT32_ID, &
    SIDRE_UINT64_ID     => CONDUIT_UINT64_ID, &
    SIDRE_FLOAT32_ID    => CONDUIT_FLOAT32_ID, &
    SIDRE_FLOAT64_ID    => CONDUIT_FLOAT64_ID, &
    SIDRE_CHAR8_STR_ID  => CONDUIT_CHAR8_STR_ID, &
    SIDRE_INT_ID        => CONDUIT_INT_ID, &
    SIDRE_UINT_ID       => CONDUIT_UINT_ID, &
    SIDRE_LONG_ID       => CONDUIT_LONG_ID, &
    SIDRE_ULONG_ID      => CONDUIT_ULONG_ID, &
    SIDRE_FLOAT_ID      => CONDUIT_FLOAT_ID, &
    SIDRE_DOUBLE_ID     => CONDUIT_DOUBLE_ID
use, intrinsic :: iso_c_binding, only : C_INT64_T
! splicer end module_use

! splicer begin module_top
integer, parameter :: MAXNAMESIZE = 128

integer, parameter :: SIDRE_IndexType = C_INT64_T

integer, parameter :: invalid_index = -1_SIDRE_IndexType
! splicer end module_top


# SIDRE_create_fortran_allocatable_view is not in api.yaml since it is not in src/core and
# only required for the fortran interface.

! splicer begin additional_interfaces
function SIDRE_create_array_view(group, name, lname, addr, type, rank, extents) &
      result(rv) bind(C,name="SIDRE_create_array_view")
    use iso_c_binding
    import SIDRE_IndexType
    type(C_PTR), value, intent(IN)     :: group
    character(kind=C_CHAR), intent(IN) :: name(*)
    integer(C_INT), value, intent(IN)  :: lname
    type(C_PTR), value,     intent(IN) :: addr
    integer(C_INT), value, intent(IN)  :: type
    integer(C_INT), value, intent(IN)  :: rank
    integer(SIDRE_IndexType), intent(IN)  :: extents(*)
    type(C_PTR) rv
end function SIDRE_create_array_view
! splicer end additional_interfaces
