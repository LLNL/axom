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

!-------------------------------------------------------------------------------
! Simple Fortran example of the quest_inout interface
!
! NOTE: assumes a "sphere.stl" file exists in the working directory.
!-------------------------------------------------------------------------------

module quest_interface_ex

#include "axom/config.hpp"

 use iso_c_binding, only : C_DOUBLE

#if defined AXOM_USE_MPI && ! defined AXOM_USE_MPIF_HEADER
  use mpi
#endif
  use axom_quest
  implicit none

#if defined AXOM_USE_MPI && defined AXOM_USE_MPIF_HEADER
  include 'mpif.h'
#endif

  integer :: num_pts = 10, rank, nprocs, ierr
  character(LEN=:), allocatable :: FileName
  integer :: dim = 3

contains

  !-------------------------------------------------------------------------------
  ! random_double( lb, ub )
  !
  ! Returns a random double, d, between lowerbound lb and upperbound ub
  !-------------------------------------------------------------------------------
  function random_double(lb,ub ) result(d)
    implicit none

    real(C_DOUBLE), value, intent(in) ::  lb, ub
    real(C_DOUBLE) :: d

    call random_number(d)
    d =  (ub-lb) * d + lb

  end function random_double

  subroutine init_mpi()
    !...initialize MPI
#ifdef AXOM_USE_MPI
    call mpi_init( ierr )
    call mpi_comm_rank( MPI_COMM_WORLD, rank, ierr )
    call mpi_comm_size( MPI_COMM_WORLD, nprocs, ierr )
#endif
  end subroutine init_mpi

  subroutine finalize_mpi()
    !... finalize MPI
#ifdef AXOM_USE_MPI
    call mpi_finalize( ierr )
#endif
  end subroutine finalize_mpi

  !-------------------------------------------------------------------------------
  ! run_distance_queries( rank )
  !
  ! Exercises the distance query functionality of quest on a given mpi rank
  !-------------------------------------------------------------------------------
  subroutine run_distance_queries()
    implicit none

    integer :: i, ins, rc
    real(C_DOUBLE)  :: x, y, z, phi
    integer :: max_occupancy = 25
    integer :: max_levels = 10
    real(C_DOUBLE) :: lo(3), hi(3)
    character(LEN=40) :: Format = "(I2, F6.2, F6.2, F6.2, F6.2, I3)"

    if (rank == 0) then
      write(*,*) "** Running distance queries using quest fortran interface"
      write(*,*) "[rank, (x,y,z), dist, within]"
    end if


    !...initialize quest
    call quest_signed_distance_set_max_levels( max_levels )
    call quest_signed_distance_set_max_occupancy( max_occupancy )
#ifdef AXOM_USE_MPI
    rc = quest_signed_distance_init( FileName, MPI_COMM_WORLD )
#else
    rc = quest_signed_distance_init( FileName )
#endif

    if ( rc /= 0 ) then
      stop 'signed_distance_init() failed!'
    endif

    call quest_signed_distance_get_mesh_bounds( lo, hi)

    !... Create random pts and query surface distance and containment
    do i=1,num_pts + 1

      x = random_double( lo(1), hi(1) )
      y = random_double( lo(2), hi(2) )
      z = random_double( lo(3), hi(3) )

      phi = quest_signed_distance_evaluate( x, y, z )
      ins = phi / abs(phi)

      write(*,Format) rank, x, y, z, phi, ins
    end do

    !...finalize quest
    call quest_signed_distance_finalize()

  end subroutine run_distance_queries


  !-------------------------------------------------------------------------------
  ! run_containment_queries( rank )
  !
  ! Exercises the containment query functionality of quest on a given mpi rank
  !-------------------------------------------------------------------------------
  subroutine run_containment_queries()
    implicit none

    integer :: i, ins
    real(C_DOUBLE)  :: x, y, z
    integer :: unused = -1
    real(C_DOUBLE) :: mesh_bb_min(3), mesh_bb_max(3), mesh_center_of_mass(3)
    character(LEN=40) :: Format = "(I2, F6.2, F6.2, F6.2, I3)"

    if (rank == 0) then
      write(*,*) "** Running containment queries using quest fortran interface"
      write(*,*) "[rank, (x,y,z), within]"
    end if


    !...initialize quest
#ifdef AXOM_USE_MPI
    call quest_initialize( MPI_COMM_WORLD, FileName, .false., dim, unused, unused )
#else
    call quest_initialize( FileName, .false., dim, unused, unused )
#endif

    call quest_mesh_min_bounds(mesh_bb_min)
    call quest_mesh_max_bounds(mesh_bb_max)
    call quest_mesh_center_of_mass(mesh_center_of_mass)

    !... Create random pts and query surface distance and containment
    do i=1,num_pts + 1
      if (i == 1) then
        x = mesh_center_of_mass(1)
        y = mesh_center_of_mass(2)
        z = mesh_center_of_mass(3)
      else
        x = random_double( mesh_bb_min(1), mesh_bb_max(1) )
        y = random_double( mesh_bb_min(2), mesh_bb_max(2) )
        z = random_double( mesh_bb_min(3), mesh_bb_max(3) )
      end if

      ins = quest_inside(x,y,z)

      write(*,Format) rank, x, y, z, ins
    end do

    !...finalize quest
    call quest_finalize()
  end subroutine run_containment_queries

end module quest_interface_ex


!-------------------------------------------------------------------------------
! Main Program
!-------------------------------------------------------------------------------
program quest_fortran

  use quest_interface_ex

  implicit none
  integer len_name

  call init_mpi()

  call get_command_argument(1, length=len_name)
  allocate(character(len_name):: FileName)
  call get_command_argument(1, FileName)

  !...exercise the distance functions in quest module
  call run_distance_queries()

  !...exercise the containment functions in quest module
  call run_containment_queries()

  deallocate(FileName)
  call finalize_mpi()

end program
