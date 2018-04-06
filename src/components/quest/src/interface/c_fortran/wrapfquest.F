! wrapfquest.F
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
!! \file wrapfquest.F
!! \brief Shroud generated wrapper for QUEST library
!<
! splicer begin file_top
#include <axom/config.hpp>
! splicer end file_top
module axom_quest
    ! splicer begin module_use
    ! splicer end module_use
    implicit none

    ! splicer begin module_top
    ! splicer end module_top

    interface

#ifdef AXOM_USE_MPI
        subroutine c_initialize_mpi(comm, fileName, requiresDistance, &
                ndims, maxElements, maxLevels) &
                bind(C, name="QUEST_initialize_mpi")
            use iso_c_binding, only : C_BOOL, C_CHAR, C_INT
            implicit none
            integer(C_INT), value, intent(IN) :: comm
            character(kind=C_CHAR), intent(IN) :: fileName(*)
            logical(C_BOOL), value, intent(IN) :: requiresDistance
            integer(C_INT), value, intent(IN) :: ndims
            integer(C_INT), value, intent(IN) :: maxElements
            integer(C_INT), value, intent(IN) :: maxLevels
        end subroutine c_initialize_mpi
#endif

#ifdef AXOM_USE_MPI
        subroutine c_initialize_mpi_bufferify(comm, fileName, LfileName, &
                requiresDistance, ndims, maxElements, maxLevels) &
                bind(C, name="QUEST_initialize_mpi_bufferify")
            use iso_c_binding, only : C_BOOL, C_CHAR, C_INT
            implicit none
            integer(C_INT), value, intent(IN) :: comm
            character(kind=C_CHAR), intent(IN) :: fileName(*)
            integer(C_INT), value, intent(IN) :: LfileName
            logical(C_BOOL), value, intent(IN) :: requiresDistance
            integer(C_INT), value, intent(IN) :: ndims
            integer(C_INT), value, intent(IN) :: maxElements
            integer(C_INT), value, intent(IN) :: maxLevels
        end subroutine c_initialize_mpi_bufferify
#endif

#ifndef AXOM_USE_MPI
        subroutine c_initialize_serial(fileName, requiresDistance, &
                ndims, maxElements, maxLevels) &
                bind(C, name="QUEST_initialize_serial")
            use iso_c_binding, only : C_BOOL, C_CHAR, C_INT
            implicit none
            character(kind=C_CHAR), intent(IN) :: fileName(*)
            logical(C_BOOL), value, intent(IN) :: requiresDistance
            integer(C_INT), value, intent(IN) :: ndims
            integer(C_INT), value, intent(IN) :: maxElements
            integer(C_INT), value, intent(IN) :: maxLevels
        end subroutine c_initialize_serial
#endif

#ifndef AXOM_USE_MPI
        subroutine c_initialize_serial_bufferify(fileName, LfileName, &
                requiresDistance, ndims, maxElements, maxLevels) &
                bind(C, name="QUEST_initialize_serial_bufferify")
            use iso_c_binding, only : C_BOOL, C_CHAR, C_INT
            implicit none
            character(kind=C_CHAR), intent(IN) :: fileName(*)
            integer(C_INT), value, intent(IN) :: LfileName
            logical(C_BOOL), value, intent(IN) :: requiresDistance
            integer(C_INT), value, intent(IN) :: ndims
            integer(C_INT), value, intent(IN) :: maxElements
            integer(C_INT), value, intent(IN) :: maxLevels
        end subroutine c_initialize_serial_bufferify
#endif

        function c_distance_0(x, y) &
                result(SHT_rv) &
                bind(C, name="QUEST_distance_0")
            use iso_c_binding, only : C_DOUBLE
            implicit none
            real(C_DOUBLE), value, intent(IN) :: x
            real(C_DOUBLE), value, intent(IN) :: y
            real(C_DOUBLE) :: SHT_rv
        end function c_distance_0

        function c_distance_1(x, y, z) &
                result(SHT_rv) &
                bind(C, name="QUEST_distance_1")
            use iso_c_binding, only : C_DOUBLE
            implicit none
            real(C_DOUBLE), value, intent(IN) :: x
            real(C_DOUBLE), value, intent(IN) :: y
            real(C_DOUBLE), value, intent(IN) :: z
            real(C_DOUBLE) :: SHT_rv
        end function c_distance_1

        function c_inside_0(x, y) &
                result(SHT_rv) &
                bind(C, name="QUEST_inside_0")
            use iso_c_binding, only : C_DOUBLE, C_INT
            implicit none
            real(C_DOUBLE), value, intent(IN) :: x
            real(C_DOUBLE), value, intent(IN) :: y
            integer(C_INT) :: SHT_rv
        end function c_inside_0

        function c_inside_1(x, y, z) &
                result(SHT_rv) &
                bind(C, name="QUEST_inside_1")
            use iso_c_binding, only : C_DOUBLE, C_INT
            implicit none
            real(C_DOUBLE), value, intent(IN) :: x
            real(C_DOUBLE), value, intent(IN) :: y
            real(C_DOUBLE), value, intent(IN) :: z
            integer(C_INT) :: SHT_rv
        end function c_inside_1

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

        subroutine quest_finalize() &
                bind(C, name="QUEST_finalize")
            implicit none
        end subroutine quest_finalize

        ! splicer begin additional_interfaces
        ! splicer end additional_interfaces
    end interface

    interface quest_distance
        module procedure quest_distance_0
        module procedure quest_distance_1
    end interface quest_distance

    interface quest_initialize
#ifdef AXOM_USE_MPI
        module procedure quest_initialize_mpi
#endif
#ifndef AXOM_USE_MPI
        module procedure quest_initialize_serial
#endif
    end interface quest_initialize

    interface quest_inside
        module procedure quest_inside_0
        module procedure quest_inside_1
    end interface quest_inside

contains

#ifdef AXOM_USE_MPI
    subroutine quest_initialize_mpi(comm, fileName, requiresDistance, &
            ndims, maxElements, maxLevels)
        use iso_c_binding, only : C_BOOL, C_INT
        integer, value, intent(IN) :: comm
        character(*), intent(IN) :: fileName
        logical, value, intent(IN) :: requiresDistance
        logical(C_BOOL) SH_requiresDistance
        integer(C_INT), value, intent(IN) :: ndims
        integer(C_INT), value, intent(IN) :: maxElements
        integer(C_INT), value, intent(IN) :: maxLevels
        SH_requiresDistance = requiresDistance  ! coerce to C_BOOL
        ! splicer begin function.initialize_mpi
        call c_initialize_mpi_bufferify(comm, fileName, &
            len_trim(fileName, kind=C_INT), SH_requiresDistance, ndims, &
            maxElements, maxLevels)
        ! splicer end function.initialize_mpi
    end subroutine quest_initialize_mpi
#endif

#ifndef AXOM_USE_MPI
    subroutine quest_initialize_serial(fileName, requiresDistance, &
            ndims, maxElements, maxLevels)
        use iso_c_binding, only : C_BOOL, C_INT
        character(*), intent(IN) :: fileName
        logical, value, intent(IN) :: requiresDistance
        logical(C_BOOL) SH_requiresDistance
        integer(C_INT), value, intent(IN) :: ndims
        integer(C_INT), value, intent(IN) :: maxElements
        integer(C_INT), value, intent(IN) :: maxLevels
        SH_requiresDistance = requiresDistance  ! coerce to C_BOOL
        ! splicer begin function.initialize_serial
        call c_initialize_serial_bufferify(fileName, &
            len_trim(fileName, kind=C_INT), SH_requiresDistance, ndims, &
            maxElements, maxLevels)
        ! splicer end function.initialize_serial
    end subroutine quest_initialize_serial
#endif

    function quest_distance_0(x, y) &
            result(SHT_rv)
        use iso_c_binding, only : C_DOUBLE
        real(C_DOUBLE), value, intent(IN) :: x
        real(C_DOUBLE), value, intent(IN) :: y
        real(C_DOUBLE) :: SHT_rv
        ! splicer begin function.distance_0
        SHT_rv = c_distance_0(x, y)
        ! splicer end function.distance_0
    end function quest_distance_0

    function quest_distance_1(x, y, z) &
            result(SHT_rv)
        use iso_c_binding, only : C_DOUBLE
        real(C_DOUBLE), value, intent(IN) :: x
        real(C_DOUBLE), value, intent(IN) :: y
        real(C_DOUBLE), value, intent(IN) :: z
        real(C_DOUBLE) :: SHT_rv
        ! splicer begin function.distance_1
        SHT_rv = c_distance_1(x, y, z)
        ! splicer end function.distance_1
    end function quest_distance_1

    function quest_inside_0(x, y) &
            result(SHT_rv)
        use iso_c_binding, only : C_DOUBLE, C_INT
        real(C_DOUBLE), value, intent(IN) :: x
        real(C_DOUBLE), value, intent(IN) :: y
        integer(C_INT) :: SHT_rv
        ! splicer begin function.inside_0
        SHT_rv = c_inside_0(x, y)
        ! splicer end function.inside_0
    end function quest_inside_0

    function quest_inside_1(x, y, z) &
            result(SHT_rv)
        use iso_c_binding, only : C_DOUBLE, C_INT
        real(C_DOUBLE), value, intent(IN) :: x
        real(C_DOUBLE), value, intent(IN) :: y
        real(C_DOUBLE), value, intent(IN) :: z
        integer(C_INT) :: SHT_rv
        ! splicer begin function.inside_1
        SHT_rv = c_inside_1(x, y, z)
        ! splicer end function.inside_1
    end function quest_inside_1

    ! splicer begin additional_functions
    ! splicer end additional_functions

end module axom_quest
