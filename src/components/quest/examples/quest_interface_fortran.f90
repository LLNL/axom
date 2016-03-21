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
! Simple Fortran example of using quest
!
! NOTE: assumes a "sphere.stl" file exists in the working directory.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! random_double( i )
!
! Returns a random double, d, given by d = i*0.5
!-------------------------------------------------------------------------------
function random_double( i ) result(d)
  implicit none
  integer, intent(in) :: i
  real*8 :: d
  d = i*0.5
end function random_double

!-------------------------------------------------------------------------------
! Main Program
!-------------------------------------------------------------------------------
program quest_fortran

  use mpi
  use quest_mod

  implicit none

  integer :: i, ierr, nprocs, rank
  real*8  :: num, x, y, z, phi, random_double

  !...initialize MPI
  call mpi_init( ierr )
  call mpi_comm_rank( MPI_COMM_WORLD, rank, ierr )
  call mpi_comm_size( MPI_COMM_WORLD, nprocs, ierr )

  !...initialize quest
  call quest_initialize( MPI_COMM_WORLD, 'sphere.stl', 3, 25, 10 )

  !... loop 10 times and get a random x,y,z to call quest_distance(x,y,z)
  do i=1,10

    x = random_double( i )
    y = random_double( i )
    z = random_double( i )

    phi = quest_distance( x, y, z )
    write(*,*) rank, x, y, z, phi

  end do

  !...finalize quest
  call quest_finalize()

  !... finalize MPI
  call mpi_finalize( ierr )

end program
