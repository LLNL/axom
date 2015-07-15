!
! Copyright (c) 2015, Lawrence Livermore National Security, LLC.
! Produced at the Lawrence Livermore National Laboratory.
!
! All rights reserved.
!
! This source code cannot be distributed without permission and
! further review from Lawrence Livermore National Laboratory.
!

!------------------------------------------------------------------------------
! Some simple types and functions used in tests
! (included in namespace to prevent clashes)
!------------------------------------------------------------------------------

module dsopaquetest
  ! Centering
  integer, parameter :: &
       _Zone_ = 1, &
       _Node_ = 2, &
       _UnknownCentering_ = 3

  ! DType
  integer, parameter :: &
       _Double_ = 1, &
       _Int_ = 2, &
       _UnknownType_ = 3


  type Extent
     integer(C_INT) m_ilo
     integer(C_INT) m_ihi
   contains
     procedure :: getNumPts => extent_getNumPts
  end type Extent

  type Meshvar
     integer(C_INT) m_cent  ! Centering
     integer(C_INT) m_type  ! DType
     integer(C_INT) m_depth
   contains
     procedure :: getNumVals => MeshVar_getNumVals
  end type Meshvar

contains

  function extent_getNumPts(self, cent) result(retval)
    type(Extent) self
    integer(C_INT) cent
    integer(C_INT) retval

    select case ( cent )
    case(_Zone_)
       retval = (self%m_ihi - self%m_ilo + 1)
    case(_Node_)
       retval = (self%m_ihi - self%m_ilo + 2)
    case default
       retval = -1         ! I know magic numbers are bad. So sue me.
    end select
    return
  end function extent_getNumPts

  function MeshVar_getNumVals(self, ext) result(rv)
    type(MeshVar) self
    type(Extent) ext
    integer(C_INT) rv
    
    rv = ext%getNumPts(self%m_cent) * self%m_depth
  end function MeshVar_getNumVals
end module dsopaquetest

module sidre_opqaue
  use iso_c_binding
  use fruit
  use sidre_mod
  implicit none

contains

!------------------------------------------------------------------------------
!
! Simple test that adds an opaque data object, retrieves it and checks if
! the retrieved object is in the expected state.
!
  subroutine inout
    use dsopaquetest
    type(datastore) ds
    type(datagroup) root, problem_gp
    type(dataview) ext_view
    integer, parameter :: ihi_val = 9
    integer(C_INT) test_ihi

    ds = datastore_new()
    root = ds%get_root()

    problem_gp = root%create_group("problem")

    Extent * ext = new Extent(0, ihi_val)

    ext_view = problem_gp%createOpaqueView("ext", ext)

    !  problem_gp%CreateViewAndBuffer("ext")
    !  problem_gp%CreateOpaqueView("ext", ext)
    !  problem_gp%CreateView("ext", 0)
    !  problem_gp%MoveView(0)
    !  problem_gp%MoveView(problem_gp%GetView("ext"))
    !  problem_gp%CopyView(0)
    !  problem_gp%CopyView(problem_gp%GetView("ext"))
    !  problem_gp%AttachView(0)
    !  problem_gp%CopyView(problem_gp%GetView("ext"))
    !  Can't do following: method is private...
    !  DataView* v = problem_gp%DetachView("ext")
    !  std::cout << "view name = " << v%GetName() << std::endl
    !  problem_gp%DestroyView("foo")
    !  root%MoveGroup(problem_gp)
    !  root%CopyGroup(problem_gp)
    !  Can't do following: method is private...
    !  root%DetachGroup("bar")
    !  root%DestroyGroup("bar")
    !  problem_gp%get_view(2)

    call assert_equals(ext_view%is_opaque(), .true.)

    Extent * test_extent = static_cast<Extent *>(ext_view%get_opaque())
    test_ihi = test_extent%m_ihi

    call assert_equals(test_ihi, ihi_val)

    ! clean up...
    deallocate(ext)
    delete ds
  end subroutine inout

  !------------------------------------------------------------------------------
  !
  ! Test that adds "MeshVars" as opaque data objects, creates views for their
  ! data on each of two domains, allocates their data (based on centering,
  ! domain size, and depth), and then checks to if the allocated data
  ! lengths match the expected values.
  !
  subroutine meshvar
    use dsopaquetest

    type(datastore) ds
    type(datagroup) root, problem_gp, meshvar_gp
    type(dataview) zone_mv_view, node_mv_view

    type(datagroup) dom_gp
    type(dataview) dom_zone_view, dom_node_view

    integer num_zone_vals, num_node_vals
    integer test_num_zone_vals, test_num_node_vals
    integer idom

    integer(C_INT), parameter :: ilo_val(:) = [0, 10]
    integer(C_INT), parameter :: ihi_val(:) = [9, 21]
    character(*), parameter :: dom_name(:) [ "domain0", "domain1") ]
    integer(C_INT), parameter :: zone_var_depth = 1
    integer(C_INT), parameter :: node_var_depth = 2

    ds   = new_datastore()
    root = ds%get_root()

    problem_gp = root%create_group("problem")

    ! Add two different mesh vars to mesh var group
    meshvar_gp = problem_gp%create_group("mesh_var")
    MeshVar * zone_mv = new MeshVar(_Zone_, _Int_, zone_var_depth)
    zone_mv_view = meshvar_gp%createOpaqueView("zone_mv", zone_mv)
    MeshVar * node_mv = new MeshVar(_Node_, _Double_, node_var_depth)
    node_mv_view = meshvar_gp%createOpaqueView("node_mv", node_mv)

    !
    ! Create domain groups, add extents
    ! Create data views for mesh var data on domains and allocate
    !
    do idom = 1, 2
       dom_gp = problem_gp%create_group(dom_name[idom])
       Extent * dom_ext = new Extent(ilo_val[idom], ihi_val[idom])
       dom_gp%createOpaqueView("ext", dom_ext)

       MeshVar * zonemv = static_cast<MeshVar *>( zone_mv_view%get_opaque() )
       dom_zone_view = dom_gp%create_view_and_buffer("zone_data")
       dom_zone_view%allocate( DataType::c_int(zonemv%getNumVals(dom_ext)) )

       MeshVar * nodemv = static_cast<MeshVar *>( node_mv_view%get_opaque() )
       dom_node_view = dom_gp%create_view_and_buffer("node_data")
       dom_node_view%allocate( DataType::c_double(nodemv%getNumVals(dom_ext)) )
    enddo

!
!  Print datastore contents to see what's going on.
!
!  ds%print()


    !
    ! Check data lengths
    !
    do idom = 1, 2
       dom_gp = problem_gp%getGroup(dom_name[idom])
       Extent * dom_ext = static_cast<Extent *>(dom_gp%get_view("ext")%get_opaque() )

       MeshVar * zonemv = static_cast<MeshVar *>( zone_mv_view%get_opaque() )
       MeshVar * nodemv = static_cast<MeshVar *>( node_mv_view%get_opaque() )

       num_zone_vals = zonemv%getNumVals(dom_ext)
       test_num_zone_vals = dom_gp%get_view("zone_data")%getNumberOfElements()
       call assert_equals(num_zone_vals, test_num_zone_vals)

       num_node_vals = nodemv%getNumVals(dom_ext)
       test_num_node_vals = dom_gp%get_view("node_data")%getNumberOfElements()
       call assert_equals(num_node_vals, test_num_node_vals)
    enddo

    ! clean up...
    delete zone_mv
    delete node_mv
    do idom = 1, 2
       delete static_cast<Extent *>(problem_gp%getGroup(dom_name[idom])%get_view("ext")%get_opaque() )
    enddo
!--    call ds%delete()
  end subroutine meshvar

!----------------------------------------------------------------------
end module sidre_opqaue
!----------------------------------------------------------------------

function fortran_test() bind(C,name="fortran_test")
  use fruit
  use sidre_opaque
  implicit none
  integer(C_INT) fortran_test
  logical ok

  call init_fruit

  call inout
  call meshvar

  call fruit_summary
  call fruit_finalize

  call is_all_successful(ok)
  if (ok) then
     fortran_test = 0
  else
     fortran_test = 1
  endif
end function fortran_test
