!
! copyright (c) 2015, lawrence livermore national security, llc.
! produced at the lawrence livermore national laboratory.
!
! all rights reserved.
!
! this source code cannot be distributed without permission and
! further review from lawrence livermore national laboratory.
!

program spio_basis_write_read

  use iso_c_binding
  use sidre_mod
  use spio_mod
  implicit none

  include 'mpif.h'

  integer mpierr
  integer num_files
  integer return_val

  type(datastore) ds, ds2
  type(datagroup) root, root2
  type(datagroup) flds, flds2
  type(datagroup) ga, gb
  type(dataview)  view1, view2

!  type(io_manager) writer, reader

  call mpi_init(mpierr)

  ds = datastore_new()
  root = ds%get_root()

  flds = root%create_group("fields")
  flds2 = root%create_group("fields2")

  ga = flds%create_group("a")
  gb = flds2%create_group("b")

  view1 = ga%create_view_scalar_int("i0", 101)
  view2 = gb%create_view_scalar_int("i1", 404)

  num_files = 1
!  write = writer_new(mpi_comm_world)

!  call writer%write(root, num_files, "F_out_spio_basic_write_read", "conduit_hdf5")

  ds2 = datastore_new()

!  reader = reader_new(mpi_comm_world)

!  root2 = ds2%get_root()
!  call reader%read(root2, "F_out_spio_basic_write_read.root")

!  return_val = 0
!  if (roo2%is_equivalent_to(root) == .false.) then
!     return_val = 1 
!  endif

!  int testvalue =
!    call ds%get_root()%get_group("fields")%get_group("a")%get_view("i0")%get_data()
!  int testvalue2 =
!    call ds2%get_root()%get_group("fields")%get_group("a")%get_view("i0")%get_data()
!
!  if (testvalue != testvalue2) {
!    return_val = 1
!  }
!
!  testvalue =
!    call ds%get_root()%get_group("fields2")%get_group("b")%get_view("i1")%get_data()
!  testvalue2 =
!    call ds2%get_root()%get_group("fields2")%get_group("b")%get_view("i1")%get_data()
!
!  if (testvalue != testvalue2) {
!    return_val = 1
!  }

  call ds%delete()
  call ds2%delete()

  call mpi_finalize(mpierr)

  call exit(return_val)
end program spio_basis_write_read
