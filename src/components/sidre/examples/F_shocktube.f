! gfortran -g -Wall -fbounds-check
!*************************************************************************
! Program:  ShockTube.C
! Purpose:  1D shock tube, split flux Euler equations
!
!         | m  |            |    mv    |
!     Q = | mv |        F = | mv^2 + P |
!         | E  |            |  v(E+P)  |
!
!     P = (gamma - 1.0)[E - 0.5 mv^2 ]
!
!             Cp
!     gamma = --     m = mass/volume   v = velocity
!             Cv
!
!     All quantiites are non-dimensionalized.
!
!     @Q   @F    @Q   @F @Q
!   -- + -- =  -- + -- -- = 0
!     @t   @x    @t   @Q @x
!
!*************************************************************************

!*************************************************************************
! Subroutine:  main
! Purpose   :  Simulate a 1D Shock Tube using split flux Euler formulation
!*************************************************************************

program main

  use iso_c_binding
  use sidre_mod
  implicit none

  integer, parameter :: IUPWIND = 1, IDOWNWIND = 2
  real(8), parameter :: gammaa = sqrt(2.0d0)  ! M_SQRT2
  real(8), parameter :: gammaaInverse = 1.0d0 /gammaa ! M_SQRT1_2


 !--// extern void DumpUltra(DataGroup * const prob)

  ! We should be able to parallelize pretty easily by
  ! adding an MPI_Init() here, modifying the setup slightly,
  ! adding a communication subroutine in the main loop,
  ! and calling MPI_Finalize() at the end of main()

  type(datastore) ds
  type(datagroup) root, prob

  integer currCycle
  integer numCyclesPerDump
  integer numUltraDumps
  integer numTotalCycles
  integer numElems
  integer numFaces
  integer numTubeElems
  real(8) time
  integer cycle
  real(8) dx, dt
  real(8) pressureRatio
  real(8) densityRatio
!  real(8), allocatable :: mass(:), momentum(:), energy(:), pressure(:)
!  integer, allocatable :: elemToFace(:,:), mapToElems(:)
!  integer, allocatable :: faceToElem(:,:)
!  real(8), allocatable :: F0(:), F1(:), F2(:)

  ds   = datastore_new()
  root = ds%get_root()
  prob = root%create_group("problem")

  call GetUserInput(prob)
  call CreateShockTubeMesh(prob)
  call InitializeShockTube(prob)


  ! use a reference when you want to update the param directly
!  currCycle = integer* const currCycle = prob->GetDataObject("cycle")->GetData<int*>()
!  numTotalCycles = *(prob->GetDataObject("numTotalCycles")->GetData<int*>())
!  dumpInterval = *(prob->GetDataObject("numCyclesPerDump")->GetData<int*>())

!!!  for (*currCycle = 0 *currCycle < numTotalCycles ++(*currCycle) )
  do currCycle = 1, numTotalcycles
    ! dump the ultra file, based on the user chosen attribute mask
    if ( mod(currCycle, numCyclesperDump) == 0) then
      call DumpUltra(prob)
   endif

    call ComputeFaceInfo(prob)
    call UpdateElemInfo(prob)
  enddo

  call DumpUltra(prob) ! One last dump

  call datastore_delete(ds)

contains


  subroutine CreateScalarIntBufferViewAndSetVal(grp, name, value)
    type(datagroup) grp
    character(*) name
    integer(C_INT) value
    type(dataview) tmpview

    tmpview = grp%create_view_and_buffer(name, ATK_C_INT_T, 1)
    call tmpview%allocate()
    call tmpview%set_value(value)
  end subroutine CreateScalarIntBufferViewAndSetVal


  subroutine CreateScalarFloatBufferViewAndSetVal(grp, name, value)
    type(datagroup) grp
    character(*) name
    real(C_DOUBLE) value
    type(dataview) tmpview

    tmpview = grp%create_view_and_buffer(name, ATK_C_DOUBLE_T, 1)
    call tmpview%allocate()
    call tmpview%set_value(value)
  end subroutine CreateScalarFloatBufferViewAndSetVal

!*************************************************************************
! Subroutine:  GetUserInput
! Purpose   :  Ask for control and output information
!*************************************************************************

subroutine GetUserInput(prob)

  type(datagroup), intent(IN) :: prob
  real(8) :: pratio, dratio

  !********************************!
  ! Get mesh info, and create mesh !
  !********************************!
  print *, "How many zones for the 1D shock tube? "
  !--//read *, numElems
  numElems = 10

  numElems = numElems + 2 ! add an inflow and outflow zone
  numFaces = numElems - 1

  call CreateScalarIntBufferViewAndSetVal(prob, "numElems", numElems)

  call CreateScalarIntBufferViewAndSetVal(prob, "numFaces", numFaces)

  !******************!
  ! Get physics info !
  !******************!
  pratio = -1.0
  dratio = -1.0

!--#if 0
!--  print *, "What cfl number would you like to use? "
!--  read *, cfl
!-#endif

  do while (pratio < 0.0 .or. pratio > 1.0)
     print *, "What pressure ratio would you like (0 <= x <= 1)? "
     !--//read *, pratio)
     pratio = 0.5
  enddo

  do while (dratio < 0.0 .or. dratio > 1.0)
     print *, "What density ratio would you like (0 <= x <= 1)? "
     !--//read *, dratio
     dratio = 0.5
  enddo


  call CreateScalarFloatBufferViewAndSetVal(prob, "pressureRatio", pratio)

  call CreateScalarFloatBufferViewAndSetVal(prob, "densityRatio", dratio)
!!  pressureRatio = pratio
!!  densityRatio = dratio

  !******************!
  ! Get output  info !
  !******************!
  print *, "How many Ultra dumps would you like? "
!--//    read *, numUltraDumps
  numUltraDumps = 10

  print *, "How many cycles per Ultra dump would you like? "
!--//    read *, numCyclesPerDump
  numCyclesPerDump = 10

  call CreateScalarIntBufferViewAndSetVal(prob, "numUltraDumps", numUltraDumps)

  call CreateScalarIntBufferViewAndSetVal(prob, "numCyclesPerDump", numCyclesPerDump)

  call CreateScalarIntBufferViewAndSetVal(prob, "numTotalCycles", numUltraDumps * numCyclesPerDump)

!!  numTotalCycles = numUltraDumps * numCyclesPerDump

  return
end subroutine GetUserInput

!**************************************************************************
! Subroutine:  CreateShockTubeMesh
! Purpose   :  Build an empty mesh for the shock tube
!
!
! Gaps between elements are faces
! |
! -------------------------------
! |   |   |             |   |   |
!
! ### ### ### ###       ### ### ### ###
! ### ### ### ###  ...  ### ### ### ###  <--- 1D Shock tube model
! ### ### ### ###       ### ### ### ###
!
! |  |                           |  |
! |  -----------------------------  |
! Inflow           |               Outflow
! Element      Tube Elements       Element
!
!
!**************************************************************************

subroutine CreateShockTubeMesh(prob)
  type(datagroup), intent(IN) :: prob
  type(datagroup) elem, face, tube
  type(datagroup) ingrp, outgrp  ! XXX unused
  type(dataview) tmpview
  type(dataview) mapToElemsView
  type(dataview) faceToElemView
  type(dataview) elemToFaceView
  integer(C_INT), pointer :: mapToElems(:)
  integer(C_INT), pointer :: faceToElem(:)
  integer(C_INT), pointer :: elemToFace(:)
  integer i, k
  integer numElems
  integer numFaces
  integer inflow(1)
  integer outflow(1)

  tmpview = prob%get_view("numElems")
  numElems = tmpview%get_value_int()
  tmpview = prob%get_view("numFaces")
  numFaces = tmpview%get_value_int()

  ! create element and face classes

  elem = prob%create_group("elem")
  face = prob%create_group("face")

  ! set up some important views

  inflow(1) = 1 ! identify inflow elements
!--//  elem->viewCreate("inflow", new IndexSet(1, inflow))
  ingrp = elem%create_group("inflow")

  outflow(1) = numElems ! identify outflow elements
!--//  elem->viewCreate("outflow", new IndexSet(1, outflow))
  outgrp = elem%create_group("outflow")

  ! identify shock tube elements - set up basic map
!--//  View *tube = elem->viewCreate("tube", new IndexSet(numElems - 2))
  ! Shift IndexSet indices to identify shocktube elements
  ! (shocktube element numbers are one through numElems-1 inclusive)
!--//  tube->indexSet()->shift(1)
  numTubeElems = numElems - 2
  tube = elem%create_group("tube")

  mapToElemsView = tube%create_view_and_buffer("mapToElems", ATK_C_INT_T, numTubeElems)
  call mapToElemsView%allocate()

  call mapToElemsView%get_value(mapToElems)
!--  allocate(mapToElem(nuMTubeElems))

  do k=1, numTubeElems
    mapToElems(k) = k + 1
  enddo

  ! Set up some important data relations

  ! Each face connects to two elements

  faceToElemView = tube%create_view_and_buffer("faceToElem", ATK_C_INT_T, 2*numFaces)
  call faceToElemView%allocate()

  call faceToElemView%get_value(faceToElem)
!--  allocate(faceToElem(2, numFaces))

  do i=1, numFaces
    faceToElem((i-1) * 2 + IUPWIND) = i
    faceToElem((i-1) * 2 + IDOWNWIND) = i + 1
  enddo

  ! Each element connects to two faces
!--//  Relation &elemToFace = *tube%relationCreate("elemToFace", 2)
!  dims(0) = numElems
  elemToFaceView = tube%create_view_and_buffer("elemToFace", ATK_C_INT_T, 2*numElems);
  call elemToFaceView%allocate()
  call elemToFaceView%get_value(elemToFace)
!!  allocate(elemToFace(2, numElems))

  do i=1, numElems
    elemToFace((i-1) * 2 + IUPWIND) = i ! same map as above by coincidence
    elemToFace((i-1) * 2 + IDOWNWIND) = i + 1
  enddo

  return
end subroutine CreateShockTubeMesh

!************************************************************************
! Subroutine:  InitializeShockTube
! Purpose   :  Populate the mesh with values
!************************************************************************

subroutine InitializeShockTube(prob)
  type(datagroup), intent(IN) :: prob
  integer i

  type(datagroup) elem, face
  type(dataview) tmpview
  integer startTube
  integer endTube
  integer midTube
  real(8) massInitial
  real(8) momentumInitial
  real(8) pressureInitial
  real(8) energyInitial
  real(8) dratio, pratio
  real(C_DOUBLE), pointer :: mass(:), momentum(:), energy(:), pressure(:)
  real(C_DOUBLE) dx

  ! These were created in GetUserInput()
  elem = prob%get_group("elem")
  face = prob%get_group("face")

  tmpview = prob%get_view("numElems")
  numElems = tmpview%get_value_int()

  tmpview = prob%get_view("numFaces")
  numFaces = tmpview%get_value_int()

  ! Create element centered quantities

  tmpview = elem%create_view_and_buffer("mass")
  call tmpview%allocate(ATK_C_DOUBLE_T, numElems)
  call tmpview%get_value(mass)

  tmpview = elem%create_view_and_buffer("momentum", ATK_C_DOUBLE_T, numElems)
  call tmpview%allocate()
  call tmpview%get_value(momentum)

  tmpview = elem%create_view_and_buffer("energy", ATK_C_DOUBLE_T, numElems)
  call tmpview%allocate()
  call tmpview%get_value(energy)

  tmpview = elem%create_view_and_buffer("pressure", ATK_C_DOUBLE_T, numElems)
  call tmpview%allocate()
  call tmpview%get_value(pressure)
!!  allocate(mass(numElems))
!!  allocate(momentum(numElems))
!!  allocate(energy(numElems))
!!  allocate(pressure(numElems))

  ! Create face centered quantities
  tmpview = elem%create_view_and_buffer("F0", ATK_C_DOUBLE_T, numFaces)
  call tmpview%allocate()
  tmpview = elem%create_view_and_buffer("F1", ATK_C_DOUBLE_T, numFaces)
  call tmpview%allocate()
  tmpview = elem%create_view_and_buffer("F2", ATK_C_DOUBLE_T, numFaces)
  call tmpview%allocate()
!!  allocate(F0(numFaces))
!!  allocate(F1(numFaces))
!!  allocate(F2(numFaces))

!--//  face%fieldCreateReal("F", 3) ! mv, mv^2+P, and v(E+P)

  ! Fill left half with high pressure, right half with low pressure
  startTube = 1
  endTube = numElems  ! elem%GetDataShape().m_dimensions(0)
  midTube = endTube / 2

  ! Non-dimensionalized reference values
  massInitial = 1.0d0
  momentumInitial = 0.0d0
  pressureInitial = gammaaInverse
  energyInitial = pressureInitial / (gammaa - 1.0d0)

  ! Initialize zonal quantities
  do i=startTube, midTube
!  for (i = startTube i < midTube ++i)
    mass(i) = massInitial
    momentum(i) = momentumInitial
    pressure(i) = pressureInitial
    energy(i) = energyInitial
  enddo

  ! adjust parameters for low pressure portion of tube
!  dratio = *(prob%GetDataObject("densityRatio")%GetData<real(8)*>())
!  pratio = *(prob%GetDataObject("pressureRatio")%GetData<real(8)*>())
  dratio = densityRatio
  pratio = pressureRatio

  massInitial = massInitial * dratio
  pressureInitial = pressureInitial * pratio
  energyInitial = pressureInitial / (gammaa - 1.0d0)

  do i=midTube, endTube
!  for (i = midTube i < endTube ++i)
    mass(i) = massInitial
    momentum(i) = momentumInitial
    pressure(i) = pressureInitial
    energy(i) = energyInitial
  enddo

  ! Create needed time info
  call CreateScalarFloatBufferViewAndSetVal(prob, "time", 0.0d0)
  call CreateScalarIntBufferViewAndSetVal(prob, "cycle", 0)

  call CreateScalarFloatBufferViewAndSetVal(prob, "dx", 1.0d0 / endTube)
  tmpview = prob%get_view("dx")
  dx = tmpview%get_value_double()
  call CreateScalarFloatBufferViewAndSetVal(prob, "dt", 0.4d0 * dx)

  time = 0.0d0
  cycle = 0
  dx = 1.0d0 / endTube
  dt = 0.4d0 * dx

  return
end subroutine InitializeShockTube

!*************************************************************************
! Subroutine:  ComputeFaceInfo
! Purpose   :  Compute F quantities at faces.
!
!  @F   @F0   @F1   @F2
!  -- = --- + --- + ---
!  @x   @x    @x    @x
!
!  Calculate F0, F1 and F2 at the face centers.
!
!*************************************************************************

subroutine ComputeFaceInfo(prob)
  type(datagroup), intent(IN) :: prob

  type(datagroup) face, elem
  type(dataview) tmpview

  integer i
  integer upWind, downWind, contributor
  real(8) massf, momentumf, energyf, pressuref, c, v, ev, cLocal

  real(C_DOUBLE), pointer :: F0(:), F1(:), F2(:)
  integer(C_INT), pointer :: faceToElem(:)
  integer(C_INT) numFaces
  real(C_DOUBLE), pointer :: mass(:), momentum(:), energy(:)

  face = prob%get_group("face")

  tmpview = face%get_view("faceToElem")
  call tmpview%get_value(faceToElem)

  tmpview = face%get_view("F0")
  call tmpview%get_value(F0)
  tmpview = face%get_view("F1")
  call tmpview%get_value(F1)
  tmpview = face%get_view("F2")
  call tmpview%get_value(F2)

  numFaces = tmpview%get_number_of_elements();

  elem = prob%get_group("elem")
  tmpview = elem%get_view("mass")
  call tmpview%get_value(mass)
  tmpview = elem%get_view("momentum")
  call tmpview%get_value(momentum)
  tmpview = elem%get_view("energy")
  call tmpview%get_value(energy)

  do i=1, numFaces
     ! each face has an upwind and downwind element.
     upWind = faceToElem((i-1) * 2 + IUPWIND)  ! upwind element
     downWind = faceToElem((i-1) * 2 + IDOWNWIND)  ! downwind element

     ! calculate face centered quantities
     massf = 0.5d0 * (mass(upWind) + mass(downWind))
     momentumf = 0.5d0 * (momentum(upWind) + momentum(downWind))
     energyf = 0.5d0 * (energy(upWind) + energy(downWind))
     pressuref = (gammaa - 1.0) * (energyf - 0.5d0 * momentumf * momentumf / massf)
     c = sqrt(gammaa * pressuref / massf)
     v = momentumf / massf

     ! Now that we have the wave speeds, we might want to
     ! look for the max wave speed here, and update dt
     ! appropriately right before leaving this function.
     ! ...

     ! OK, calculate face quantities

     F0(i) = 0.0d0
     F1(i) = 0.0d0
     F2(i) = 0.0d0

     if (v >= 0.0d0) then
        contributor = upWind
     else
        contributor = downWind
     endif
     massf = mass(contributor)
     momentumf = momentum(contributor)
     energyf = energy(contributor)
     pressuref = energyf - 0.5d0 * momentumf * momentumf / massf
     ev = v * (gammaa - 1.0d0)

     F0(i) = F0(i) + ev * massf
     F1(i) = F1(i) + ev * momentumf
     F2(i) = F2(i) + ev * (energyf - pressuref)
     
     if (v + c >= 0.0d0) then
        contributor = upWind
     else
        contributor = downWind
     endif
     massf = mass(contributor)
     momentumf = momentum(contributor)
     energyf = energy(contributor)
     pressuref = (gammaa - 1.0d0) * (energyf - 0.5d0 * momentumf * momentumf / massf)
     ev = 0.5d0 * (v + c)
     cLocal = sqrt(gammaa * pressuref / massf)

     F0(i) = F0(i) + ev * massf
     F1(i) = F1(i) + ev * (momentumf + massf * cLocal)
     F2(i) = F2(i) + ev * (energyf + pressuref + momentumf * cLocal)

     if (v - c >= 0.0d0) then
        contributor = upWind
     else
        contributor = downWind
     endif
     massf = mass(contributor)
     momentumf = momentum(contributor)
     energyf = energy(contributor)
     pressuref = (gammaa - 1.0d0) * (energyf - 0.5d0 * momentumf * momentumf / massf)
     ev = 0.5d0 * (v - c)
     cLocal = sqrt(gammaa * pressuref / massf)

     F0(i) = F0(i) + ev * massf
     F1(i) = F1(i) + ev * (momentumf - massf * cLocal)
     F2(i) = F2(i) + ev * (energyf + pressuref - momentumf * cLocal)
  enddo

  return
end subroutine ComputeFaceInfo

!*************************************************************************
! Subroutine:  UpdateElemInfo
! Purpose   :  Q(elem) = Q(elem) + deltaQ(elem)
!
!  deltaQ(elem) = - (F(downWindFace) - F(upWindFace)) * dt / dx
!
!*************************************************************************

subroutine UpdateElemInfo(prob)
  type(datagroup), intent(IN) :: prob
  type(datagroup) elem, face, tube
  type(dataview) tmpview

  real(C_DOUBLE) dx, dt, time
  integer(C_INT), pointer :: elemToFace(:)
  integer i
  integer elemIdx, upWind, downWind
  integer(C_INT), pointer :: is(:)
  real(C_DOUBLE), pointer :: F0(:), F1(:), F2(:)
  real(C_DOUBLE), pointer :: mass(:), momentum(:), energy(:), pressure(:)

  ! get the element quantities we want to update
  elem = prob%get_group("elem")

  tmpview = elem%get_view("mass")
  call tmpview%get_value(mass)
  tmpview = elem%get_view("momentum")
  call tmpview%get_value(momentum)
  tmpview = elem%get_view("energy")
  call tmpview%get_value(energy)
  tmpview = elem%get_view("pressure")
  call tmpview%get_value(pressure)

  ! focus on just the elements within the shock tube
  tube = elem%get_group("tube")

  tmpview = tube%get_view("elemToFace")
  call tmpview%get_value(elemToFace)

!--//  Relation &elemToFace = *tube%relation("elemToFace")
!  integer numTubeElems = tube%GetDataShape().m_dimensions(0)

!--//  int *is = tube%map()
  tmpview = tube%get_view("mapToElems")
  call tmpview%get_value(is)


  ! The element update is calculated as the flux between faces
  face = prob%get_group("face")
  tmpview = face%get_view("F0")
  call tmpview%get_value(F0)
  tmpview = face%get_view("F1")
  call tmpview%get_value(F1)
  tmpview = face%get_view("F2")
  call tmpview%get_value(F2)

  tmpview = prob%get_view("dx")
  dx = tmpview%get_value_double()
  tmpview = prob%get_view("dt")
  dt = tmpview%get_value_double()
  tmpview = prob%get_view("time")
  time = tmpview%get_value_double()

  do i=1, numTubeElems
     ! recalculate elements in the shocktube, don't touch inflow/outflow
     elemIdx = is(i)
     ! each element inside the tube has an upwind and downwind face
     upWind = elemToFace((i-1)*2 + IUPWIND)  ! upwind face
     downWind = elemToFace((i-1)*2 + IDOWNWIND)  ! downwind face

     mass(elemIdx) = mass(elemIdx) - (gammaaInverse * (F0(downWind) - F0(upWind)) * dt / dx)
     momentum(elemIdx) = momentum(elemIdx) - (gammaaInverse * (F1(downWind) - F1(upWind)) * dt / dx)
     energy(elemIdx) = energy(elemIdx) - (gammaaInverse * (F2(downWind) - F2(upWind)) * dt / dx)
     pressure(elemIdx) = (gammaa - 1.0d0) * (energy(elemIdx) - 0.5d0 * momentum(elemIdx) * momentum(elemIdx) / mass(elemIdx))
  enddo

  ! update the time
  time = time + dt

  return
end subroutine UpdateElemInfo

subroutine DumpUltra( prob )
  type(datagroup), intent(IN) :: prob
  character(100) fname

  type(datagroup) elem
  type(dataview) view
  integer, parameter :: fp = 8
  integer i
  integer ierr
  character(30) name
  integer length

  elem = prob%get_group("elem")

  fname = "problem"

  view = prob%get_view("cycle")

  write(fname, '(a,"_", i04, ".ult")') trim(fname), view%get_value_int()
  open(fp, file=fname, iostat=ierr)
  if (ierr /= 0) then
     print *, "Could not open file ", trim(fname), ". Aborting."
!     exit (-1)
  endif

  write(fp, '(a)') "# Ultra Plot File"
  write(fp, '(a)') "# Problem: problem"

  do i=1, prob%get_num_views()
     view = prob%get_view(i)
     length = view%get_number_of_elements()
     call view%get_name(name)
     if ( length <= 1 ) then
        select case (view%get_type_id())
        case (ATK_C_INT_T)
           write(fp, '("# ", a, " = ", i4)') name, view%get_value_int()
        case (ATK_C_DOUBLE_T)
           write(fp, '("# ", a, " = ", f12.5)') name, view%get_value_double()
        end select
    endif
  enddo


!+  for(size_t i=0  i<elem%getNumViews()  i++)
!+  {
!+    DataView * const view = elem%get_view(i)
!+    const int length = view%getSchema().dtype().number_of_elements()
!+    const std::string& name = view%get_name()
!+    write(fp, "# %s\n", name.c_str() )
!+
!+    if( view%getSchema().dtype().id() == DataType::INT32_T )
!+    {
!+      int32 const * const data = view%getValue()
!+      for ( int i=0  i<length  ++i)
!+      {
!+        write(fp, "%f %f\n", (double) i, (double) data[i])
!+      }
!+      write(fp, "\n")
!+    }
!+    else if( view%getSchema().dtype().id() == DataType::FLOAT64_T )
!+    {
!+      float64 const * const data = view%getValue()
!+      for ( int i=0  i<length  ++i)
!+      {
!+        write(fp, "%f %f\n", (double) i, (double) data[i])
!+      }
!+      write(fp, "\n")
!+    }
!+  }
!+
!+
  close(fp)

   return
 end subroutine DumpUltra

end program main
