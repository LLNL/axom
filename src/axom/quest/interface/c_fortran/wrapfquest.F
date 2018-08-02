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

#ifdef AXOM_USE_MPI
        subroutine c_signed_distance_init_mpi(file, comm) &
                bind(C, name="QUEST_signed_distance_init_mpi")
            use iso_c_binding, only : C_CHAR, C_INT
            implicit none
            character(kind=C_CHAR), intent(IN) :: file(*)
            integer(C_INT), value, intent(IN) :: comm
        end subroutine c_signed_distance_init_mpi
#endif

#ifdef AXOM_USE_MPI
        subroutine c_signed_distance_init_mpi_bufferify(file, Lfile, &
                comm) &
                bind(C, name="QUEST_signed_distance_init_mpi_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT
            implicit none
            character(kind=C_CHAR), intent(IN) :: file(*)
            integer(C_INT), value, intent(IN) :: Lfile
            integer(C_INT), value, intent(IN) :: comm
        end subroutine c_signed_distance_init_mpi_bufferify
#endif

#ifndef AXOM_USE_MPI
        subroutine c_signed_distance_init_serial(file) &
                bind(C, name="QUEST_signed_distance_init_serial")
            use iso_c_binding, only : C_CHAR
            implicit none
            character(kind=C_CHAR), intent(IN) :: file(*)
        end subroutine c_signed_distance_init_serial
#endif

#ifndef AXOM_USE_MPI
        subroutine c_signed_distance_init_serial_bufferify(file, Lfile) &
                bind(C, name="QUEST_signed_distance_init_serial_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT
            implicit none
            character(kind=C_CHAR), intent(IN) :: file(*)
            integer(C_INT), value, intent(IN) :: Lfile
        end subroutine c_signed_distance_init_serial_bufferify
#endif

        function c_signed_distance_initialized() &
                result(SHT_rv) &
                bind(C, name="QUEST_signed_distance_initialized")
            use iso_c_binding, only : C_BOOL
            implicit none
            logical(C_BOOL) :: SHT_rv
        end function c_signed_distance_initialized

        subroutine quest_signed_distance_get_mesh_bounds(lo, hi) &
                bind(C, name="QUEST_signed_distance_get_mesh_bounds")
            use iso_c_binding, only : C_DOUBLE
            implicit none
            real(C_DOUBLE), intent(OUT) :: lo(*)
            real(C_DOUBLE), intent(OUT) :: hi(*)
        end subroutine quest_signed_distance_get_mesh_bounds

        subroutine quest_signed_distance_set_dimension(dim) &
                bind(C, name="QUEST_signed_distance_set_dimension")
            use iso_c_binding, only : C_INT
            implicit none
            integer(C_INT), value, intent(IN) :: dim
        end subroutine quest_signed_distance_set_dimension

        subroutine c_signed_distance_set_closed_surface(status) &
                bind(C, name="QUEST_signed_distance_set_closed_surface")
            use iso_c_binding, only : C_BOOL
            implicit none
            logical(C_BOOL), value, intent(IN) :: status
        end subroutine c_signed_distance_set_closed_surface

        subroutine quest_signed_distance_set_max_levels(maxLevels) &
                bind(C, name="QUEST_signed_distance_set_max_levels")
            use iso_c_binding, only : C_INT
            implicit none
            integer(C_INT), value, intent(IN) :: maxLevels
        end subroutine quest_signed_distance_set_max_levels

        subroutine quest_signed_distance_set_max_occupancy(maxOccupancy) &
                bind(C, name="QUEST_signed_distance_set_max_occupancy")
            use iso_c_binding, only : C_INT
            implicit none
            integer(C_INT), value, intent(IN) :: maxOccupancy
        end subroutine quest_signed_distance_set_max_occupancy

        subroutine c_signed_distance_set_verbose(status) &
                bind(C, name="QUEST_signed_distance_set_verbose")
            use iso_c_binding, only : C_BOOL
            implicit none
            logical(C_BOOL), value, intent(IN) :: status
        end subroutine c_signed_distance_set_verbose

        function quest_signed_distance_evaluate(x, y, z) &
                result(SHT_rv) &
                bind(C, name="QUEST_signed_distance_evaluate")
            use iso_c_binding, only : C_DOUBLE
            implicit none
            real(C_DOUBLE), value, intent(IN) :: x
            real(C_DOUBLE), value, intent(IN) :: y
            real(C_DOUBLE), value, intent(IN) :: z
            real(C_DOUBLE) :: SHT_rv
        end function quest_signed_distance_evaluate

        subroutine quest_signed_distance_finalize() &
                bind(C, name="QUEST_signed_distance_finalize")
            implicit none
        end subroutine quest_signed_distance_finalize

        ! splicer begin additional_interfaces
        ! splicer end additional_interfaces
    end interface

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

    interface quest_signed_distance_init
#ifdef AXOM_USE_MPI
        module procedure quest_signed_distance_init_mpi
#endif
#ifndef AXOM_USE_MPI
        module procedure quest_signed_distance_init_serial
#endif
    end interface quest_signed_distance_init

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

#ifdef AXOM_USE_MPI
    subroutine quest_signed_distance_init_mpi(file, comm)
        use iso_c_binding, only : C_INT
        character(*), intent(IN) :: file
        integer, value, intent(IN) :: comm
        ! splicer begin function.signed_distance_init_mpi
        call c_signed_distance_init_mpi_bufferify(file, &
            len_trim(file, kind=C_INT), comm)
        ! splicer end function.signed_distance_init_mpi
    end subroutine quest_signed_distance_init_mpi
#endif

#ifndef AXOM_USE_MPI
    subroutine quest_signed_distance_init_serial(file)
        use iso_c_binding, only : C_INT
        character(*), intent(IN) :: file
        ! splicer begin function.signed_distance_init_serial
        call c_signed_distance_init_serial_bufferify(file, &
            len_trim(file, kind=C_INT))
        ! splicer end function.signed_distance_init_serial
    end subroutine quest_signed_distance_init_serial
#endif

    function quest_signed_distance_initialized() &
            result(SHT_rv)
        use iso_c_binding, only : C_BOOL
        logical :: SHT_rv
        ! splicer begin function.signed_distance_initialized
        SHT_rv = c_signed_distance_initialized()
        ! splicer end function.signed_distance_initialized
    end function quest_signed_distance_initialized

    subroutine quest_signed_distance_set_closed_surface(status)
        use iso_c_binding, only : C_BOOL
        logical, value, intent(IN) :: status
        logical(C_BOOL) SH_status
        SH_status = status  ! coerce to C_BOOL
        ! splicer begin function.signed_distance_set_closed_surface
        call c_signed_distance_set_closed_surface(SH_status)
        ! splicer end function.signed_distance_set_closed_surface
    end subroutine quest_signed_distance_set_closed_surface

    subroutine quest_signed_distance_set_verbose(status)
        use iso_c_binding, only : C_BOOL
        logical, value, intent(IN) :: status
        logical(C_BOOL) SH_status
        SH_status = status  ! coerce to C_BOOL
        ! splicer begin function.signed_distance_set_verbose
        call c_signed_distance_set_verbose(SH_status)
        ! splicer end function.signed_distance_set_verbose
    end subroutine quest_signed_distance_set_verbose

    ! splicer begin additional_functions
    ! splicer end additional_functions

end module axom_quest
