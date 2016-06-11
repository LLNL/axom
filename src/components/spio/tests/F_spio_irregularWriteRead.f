! F_spio_irregularWriteRead.f
!
! copyright (c) 2015, lawrence livermore national security, llc.
! produced at the lawrence livermore national laboratory.
!
! all rights reserved.
!
! this source code cannot be distributed without permission and
! further review from lawrence livermore national laboratory.
!

program spio_irregularWriteRead
  use sidre_mod
  use spio_mod
  implicit none

  include 'mpif.h'

  integer i, f, g
  integer mpierr
  integer num_files
  integer return_val
  integer my_rank, num_ranks, num_output
  integer num_fields, num_subgroups
  integer, pointer :: vals(:)
  character(80) name

  type(datastore) ds, ds2
  type(datagroup) root, root2
  type(datagroup) flds, flds2
  type(datagroup) sg
  type(dataview)  view

  type(iomanager) writer, reader

  call mpi_init(mpierr)

  call mpi_comm_rank(MPI_COMM_WORLD, my_rank, mpierr)
  call mpi_comm_size(MPI_COMM_WORLD, num_ranks, mpierr)

  num_output = num_ranks / 2 
  if (num_output == 0) then
     num_output = 1
  endif

! create a datastore and give it a small hierarchy of groups and views.
!
! the views are filled with repeatable nonsense data that will vary based
! on rank.
  ds = datastore_new()
  root = ds%get_root()

  num_fields = my_rank + 2

  do f = 0, num_fields-1
     write(name, "(a,i0)") "fields", f
     flds = root%create_group(name)

     num_subgroups = mod(f+my_rank,3) + 1
     do g = 0, num_subgroups-1
        write(name, "(a,i0)") "subgroup", g
        sg = flds%create_group(name)

        write(name, "(a,i0)") "view", g
        if (mod(g, 2) .ne. 0) then
           view = sg%create_view_and_allocate(name, SIDRE_INT_ID, 10+my_rank)
           call view%get_data(vals)
           do i = 0, 9+my_rank
              vals(i+1) = (i+10) * (404-my_rank-i-g-f)
           enddo
        else
           view = sg%create_view_scalar_int(name, 101*my_rank*(f+g+1))
        endif
     enddo
  enddo

  ! contents of the datastore written to files with io_manager.
  num_files = num_output
  writer = iomanager_new(MPI_COMM_WORLD)

  call writer%write(root, num_files, "F_out_spio_irregular_write_read", "conduit_hdf5")

  ! create another datastore that holds nothing but the root group.
  ds2 = datastore_new()

  ! read from the files that were written above.
  reader = iomanager_new(MPI_COMM_WORLD)

  root2 = ds2%get_root()
  call reader%read(root2, "F_out_spio_irregular_write_read.root")

  ! verify that the contents of ds2 match those written from ds.
  return_val = 0
  if (.not. root2%is_equivalent_to(root)) then
     return_val = 1
  endif

!--  for (int f = 0 f < num_fields ++f) {
!--    std::ostringstream ostream
!--    ostream << "fields" << f
!--    type(datagroup) flds
!--    flds = ds%get_root()%get_group(ostream.str())
!--    type(datagroup) flds2
!--    flds2 = ds2%get_root()%get_group(ostream.str())
!--
!--    int num_subgroups = mod(f+my_rank,3) + 1
!--    for (int g = 0 g < num_subgroups ++g) {
!--      std::ostringstream gstream
!--      gstream << "subgroup" << g
!--      type(datagroup) sg
!--      sg = flds%get_group(gstream.str())
!--      type(datagroup) sg2
!--      sg2 = flds2%get_group(gstream.str())
!--
!--      std::ostringstream vstream
!--      vstream << "view" << g
!--      if (mod(g,2) .ne. 0) then
!--
!--        type(dataview) view_orig
!--        view_orig = sg%get_view(vstream.str())
!--        type(dataview) view_restored
!--        view_restored = sg2%get_view(vstream.str())
!--
!--        int num_elems = view_orig%get_num_elements()
!--        if (view_restored%get_num_elements() != num_elems) {
!--          return_val = 1 
!--        }
!--
!--        type(int) vals_orig
!--        vals_orig = view_orig%get_data()
!--        type(int) vals_restored
!--        vals_restored = view_restored%get_data()
!--
!--        for (int i = 0 i < num_elems ++i) {
!--          if (return_val != 1) {
!--            if (vals_orig(i) != vals_restored(i)) {
!--              return_val = 1
!--            }
!--          }
!--        }
!--
!--      } else {
!--        int testvalue = sg%get_view(vstream.str())%get_data()
!--        int testvalue2 = sg2%get_view(vstream.str())%get_data()
!--
!--        if (testvalue != testvalue2) {
!--          return_val = 1
!--        }
!--      }
!--    }
!--  } 

  call ds%delete()
  call ds2%delete()

  call mpi_finalize(mpierr)

  call exit(return_val)
end program spio_irregularWriteRead

