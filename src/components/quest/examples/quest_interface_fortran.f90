!-------------------------------------------------------------------------------
!  Copyright (c) 2015, Lawrence Livermore National Security, LLC.
!  Produced at the Lawrence Livermore National Laboratory.
!
!  All rights reserved.
!
! This source code cannot be distributed without permission and
! further review from Lawrence Livermore National Laboratory.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Simple Fortran example of quest usage.
!
! NOTE: assumes a "sphere.stl" file exists in the working directory.
!-------------------------------------------------------------------------------


module quest_fortran_usage

  use mpi
  use quest_mod

  real*8 :: default_lb = -4.5
  real*8 :: default_ub = 4.5
  integer :: num_pts = 10
  character(LEN=40) :: FileName = "sphere.stl"
  integer :: dim = 3

contains

  !-------------------------------------------------------------------------------
  ! random_double( lb, ub )
  !
  ! Returns a random double, d, between lowerbound lb and upperbound ub
  !-------------------------------------------------------------------------------
  function random_double(lb,ub ) result(d)
    implicit none

    real*8, value, intent(in) ::  lb, ub
    real*8 :: d

    call random_number(d)
    d =  (ub-lb) * d + lb

  end function random_double


  !-------------------------------------------------------------------------------
  ! run_distance_queries( rank )
  !
  ! Exercises the distance query functionality of quest on a given mpi rank
  !-------------------------------------------------------------------------------
  subroutine run_distance_queries( rank )
    implicit none
    integer, intent(in) :: rank
    integer :: i, ins
    real*8  :: x, y, z, phi
    integer :: bucket_size = 25
    integer :: depth = 10
    character(LEN=40) :: Format = "(I2, F6.2, F6.2, F6.2, F6.2, I3)"

    !...initialize quest
    call quest_initialize( MPI_COMM_WORLD, FileName, .true., dim, bucket_size, depth )

    !... Create random pts and query surface distance and containment
    do i=1,num_pts
      x = random_double( default_lb, default_ub )
      y = random_double( default_lb, default_ub )
      z = random_double( default_lb, default_ub )

      phi = quest_distance( x, y, z )
      ins = quest_inside(x,y,z)

      write(*,Format) rank, x, y, z, phi, ins
    end do

    !...finalize quest
    call quest_finalize()

  end subroutine run_distance_queries


  !-------------------------------------------------------------------------------
  ! run_containment_queries( rank )
  !
  ! Exercises the containment query functionality of quest on a given mpi rank
  !-------------------------------------------------------------------------------
  subroutine run_containment_queries( rank )
    implicit none
    integer, intent(in) :: rank
    integer :: i, ins
    real*8  :: x, y, z
    integer :: unused = -1
    character(LEN=40) :: Format = "(I2, F6.2, F6.2, F6.2, I3)"

    !...initialize quest
    call quest_initialize( MPI_COMM_WORLD, FileName, .false., dim, unused, unused )

    !... Create random pts and query surface containment
    do i=1,num_pts
      x = random_double( default_lb, default_ub )
      y = random_double( default_lb, default_ub )
      z = random_double( default_lb, default_ub )

      ins = quest_inside(x,y,z)

      write(*,Format) rank, x, y, z, ins
    end do

    !...finalize quest
    call quest_finalize()
  end subroutine run_containment_queries

end module


!-------------------------------------------------------------------------------
! Main Program
!-------------------------------------------------------------------------------
program quest_fortran

  use mpi
  use quest_fortran_usage

  implicit none

  integer :: i, ierr, nprocs, rank

  !...initialize MPI
  call mpi_init( ierr )
  call mpi_comm_rank( MPI_COMM_WORLD, rank, ierr )
  call mpi_comm_size( MPI_COMM_WORLD, nprocs, ierr )

  !...exercise the distance functions in quest module
  if (rank == 0) then
    write(*,*) "** Running distance queries using quest fortran interface"
    write(*,*) "[rank, (x,y,z), dist, within]"
  end if

  call run_distance_queries( rank )


  !...exercise the containment functions in quest module
  if (rank == 0) then
    write(*,*) "** Running containment queries using quest fortran interface"
    write(*,*) "[rank, (x,y,z), within]"
  end if

  call run_containment_queries( rank )


  !... finalize MPI
  call mpi_finalize( ierr )

end program
