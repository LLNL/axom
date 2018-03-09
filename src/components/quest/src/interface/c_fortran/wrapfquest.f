! wrapfquest.f
! This is generated code, do not edit
!
! Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
!
! Produced at the Lawrence Livermore National Laboratory
!
! LLNL-CODE-741217
!
! All rights reserved.
!
! This file is part of Axom.
!
! For details about use and distribution, please read axom/LICENSE.
!
!>
!! \file wrapfquest.f
!! \brief Shroud generated wrapper for QUEST library
!<
! splicer begin file_top
! splicer end file_top
module axom_quest
    ! splicer begin module_use
    ! splicer end module_use
    implicit none

    ! splicer begin module_top
    ! splicer end module_top

    interface

        subroutine c_initialize(comm, fileName, requiresDistance, ndims, &
                maxElements, maxLevels) &
                bind(C, name="QUEST_initialize")
            use iso_c_binding, only : C_BOOL, C_CHAR, C_INT
            implicit none
            integer(C_INT), value, intent(IN) :: comm
            character(kind=C_CHAR), intent(IN) :: fileName(*)
            logical(C_BOOL), value, intent(IN) :: requiresDistance
            integer(C_INT), value, intent(IN) :: ndims
            integer(C_INT), value, intent(IN) :: maxElements
            integer(C_INT), value, intent(IN) :: maxLevels
        end subroutine c_initialize

        subroutine c_initialize_bufferify(comm, fileName, LfileName, &
                requiresDistance, ndims, maxElements, maxLevels) &
                bind(C, name="QUEST_initialize_bufferify")
            use iso_c_binding, only : C_BOOL, C_CHAR, C_INT
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
            implicit none
        end subroutine quest_finalize

        function quest_distance(x, y, z) &
                result(SHT_rv) &
                bind(C, name="QUEST_distance")
            use iso_c_binding, only : C_DOUBLE
            implicit none
            real(C_DOUBLE), value, intent(IN) :: x
            real(C_DOUBLE), value, intent(IN) :: y
            real(C_DOUBLE), value, intent(IN) :: z
            real(C_DOUBLE) :: SHT_rv
        end function quest_distance

        function quest_inside(x, y, z) &
                result(SHT_rv) &
                bind(C, name="QUEST_inside")
            use iso_c_binding, only : C_DOUBLE, C_INT
            implicit none
            real(C_DOUBLE), value, intent(IN) :: x
            real(C_DOUBLE), value, intent(IN) :: y
            real(C_DOUBLE), value, intent(IN) :: z
            integer(C_INT) :: SHT_rv
        end function quest_inside

        subroutine quest_mesh_min_bounds(coords) &
                bind(C, name="QUEST_mesh_min_bounds")
            use iso_c_binding, only : C_DOUBLE
            implicit none
            real(C_DOUBLE), intent(OUT) :: coords(*)
        end subroutine quest_mesh_min_bounds

        subroutine quest_mesh_max_bounds(coords) &
                bind(C, name="QUEST_mesh_max_bounds")
            use iso_c_binding, only : C_DOUBLE
            implicit none
            real(C_DOUBLE), intent(OUT) :: coords(*)
        end subroutine quest_mesh_max_bounds

        subroutine quest_mesh_center_of_mass(coords) &
                bind(C, name="QUEST_mesh_center_of_mass")
            use iso_c_binding, only : C_DOUBLE
            implicit none
            real(C_DOUBLE), intent(OUT) :: coords(*)
        end subroutine quest_mesh_center_of_mass

        ! splicer begin additional_interfaces
        ! splicer end additional_interfaces
    end interface

contains

    subroutine quest_initialize(comm, fileName, requiresDistance, ndims, &
            maxElements, maxLevels)
        use iso_c_binding, only : C_BOOL, C_INT
        integer, value, intent(IN) :: comm
        character(*), intent(IN) :: fileName
        logical, value, intent(IN) :: requiresDistance
        logical(C_BOOL) SH_requiresDistance
        integer(C_INT), value, intent(IN) :: ndims
        integer(C_INT), value, intent(IN) :: maxElements
        integer(C_INT), value, intent(IN) :: maxLevels
        SH_requiresDistance = requiresDistance  ! coerce to C_BOOL
        ! splicer begin function.initialize
        call c_initialize_bufferify(comm, fileName, &
            len_trim(fileName, kind=C_INT), SH_requiresDistance, ndims, &
            maxElements, maxLevels)
        ! splicer end function.initialize
    end subroutine quest_initialize

    ! splicer begin additional_functions
    ! splicer end additional_functions

end module axom_quest
