! wrapfquest.f
! This is generated code, do not edit
!
! Copyright (c) 2015, Lawrence Livermore National Security, LLC.
! Produced at the Lawrence Livermore National Laboratory.
!
! All rights reserved.
!
! This source code cannot be distributed without permission and
! further review from Lawrence Livermore National Laboratory.
!
!>
!! \file wrapfquest.f
!! \brief Shroud generated wrapper for QUEST library
!<
! splicer begin file_top
! splicer end file_top
module quest_mod
    ! splicer begin module_use
    ! splicer end module_use
    implicit none


    interface

        subroutine c_initialize(comm, fileName, requiresDistance, ndims, maxElements, maxLevels) &
                bind(C, name="QUEST_initialize")
            use iso_c_binding
            implicit none
            integer(C_INT), value, intent(IN) :: comm
            character(kind=C_CHAR), intent(IN) :: fileName(*)
            logical(C_BOOL), value, intent(IN) :: requiresDistance
            integer(C_INT), value, intent(IN) :: ndims
            integer(C_INT), value, intent(IN) :: maxElements
            integer(C_INT), value, intent(IN) :: maxLevels
        end subroutine c_initialize

        subroutine c_initialize_bufferify(comm, fileName, LfileName, requiresDistance, ndims, maxElements, maxLevels) &
                bind(C, name="QUEST_initialize_bufferify")
            use iso_c_binding
            implicit none
            integer(C_INT), value, intent(IN) :: comm
            character(kind=C_CHAR), intent(IN) :: fileName(*)
            integer(C_INT), value, intent(IN) :: LfileName
            logical(C_BOOL), value, intent(IN) :: requiresDistance
            integer(C_INT), value, intent(IN) :: ndims
            integer(C_INT), value, intent(IN) :: maxElements
            integer(C_INT), value, intent(IN) :: maxLevels
        end subroutine c_initialize_bufferify

        subroutine quest_finalize() &
                bind(C, name="QUEST_finalize")
            use iso_c_binding
            implicit none
        end subroutine quest_finalize

        function quest_distance(x, y, z) &
                result(SH_rv) &
                bind(C, name="QUEST_distance")
            use iso_c_binding
            implicit none
            real(C_DOUBLE), value, intent(IN) :: x
            real(C_DOUBLE), value, intent(IN) :: y
            real(C_DOUBLE), value, intent(IN) :: z
            real(C_DOUBLE) :: SH_rv
        end function quest_distance

        function quest_inside(x, y, z) &
                result(SH_rv) &
                bind(C, name="QUEST_inside")
            use iso_c_binding
            implicit none
            real(C_DOUBLE), value, intent(IN) :: x
            real(C_DOUBLE), value, intent(IN) :: y
            real(C_DOUBLE), value, intent(IN) :: z
            integer(C_INT) :: SH_rv
        end function quest_inside

        subroutine quest_mesh_min_bounds(coords) &
                bind(C, name="QUEST_mesh_min_bounds")
            use iso_c_binding
            implicit none
            real(C_DOUBLE), intent(OUT) :: coords(*)
        end subroutine quest_mesh_min_bounds

        subroutine quest_mesh_max_bounds(coords) &
                bind(C, name="QUEST_mesh_max_bounds")
            use iso_c_binding
            implicit none
            real(C_DOUBLE), intent(OUT) :: coords(*)
        end subroutine quest_mesh_max_bounds

        subroutine quest_mesh_center_of_mass(coords) &
                bind(C, name="QUEST_mesh_center_of_mass")
            use iso_c_binding
            implicit none
            real(C_DOUBLE), intent(OUT) :: coords(*)
        end subroutine quest_mesh_center_of_mass

        ! splicer begin additional_interfaces
        ! splicer end additional_interfaces
    end interface

contains

    subroutine quest_initialize(comm, fileName, requiresDistance, ndims, maxElements, maxLevels)
        use iso_c_binding, only : C_BOOL, C_INT
        implicit none
        integer, value, intent(IN) :: comm
        character(*), intent(IN) :: fileName
        logical, value, intent(IN) :: requiresDistance
        logical(C_BOOL) tmp_requiresDistance
        integer(C_INT), value, intent(IN) :: ndims
        integer(C_INT), value, intent(IN) :: maxElements
        integer(C_INT), value, intent(IN) :: maxLevels
        tmp_requiresDistance = requiresDistance  ! coerce to C_BOOL
        ! splicer begin initialize
        call c_initialize_bufferify(  &
            comm,  &
            fileName,  &
            len_trim(fileName, kind=C_INT),  &
            tmp_requiresDistance,  &
            ndims,  &
            maxElements,  &
            maxLevels)
        ! splicer end initialize
    end subroutine quest_initialize

    ! splicer begin additional_functions
    ! splicer end additional_functions

end module quest_mod
