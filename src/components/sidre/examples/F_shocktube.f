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
  type(datagroup) root, problem

! 'database'
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
  real(8), allocatable :: mass(:), momentum(:), energy(:), pressure(:)
  integer, allocatable :: elemToFace(:,:), mapToElems(:)
  integer, allocatable :: faceToElem(:,:)
  real(8), allocatable :: F0(:), F1(:), F2(:)

  !DataStoreNS::DataStore* const dataStore = &DATASTORE
  !DataStoreNS::DataGroup* const rootGroup = dataStore->GetRootDataGroup()
  !DataStoreNS::DataGroup* const problem = rootGroup->CreateDataGroup("problem")
  ds   = datastore_new()
  root = ds%get_root()
  problem = root%create_group("problem")

  call GetUserInput(problem)
  call CreateShockTubeMesh(problem)
  call InitializeShockTube(problem)


  ! use a reference when you want to update the param directly
!  currCycle = integer* const currCycle = problem->GetDataObject("cycle")->GetData<int*>()
!  numTotalCycles = *(problem->GetDataObject("numTotalCycles")->GetData<int*>())
!  dumpInterval = *(problem->GetDataObject("numCyclesPerDump")->GetData<int*>())

!!!  for (*currCycle = 0 *currCycle < numTotalCycles ++(*currCycle) )
  do currCycle = 1, numTotalcycles
    ! dump the ultra file, based on the user chosen attribute mask
    if ( mod(currCycle, numCyclesperDump) == 0) then
      call DumpUltra(problem)
   endif

    call ComputeFaceInfo(problem)
    call UpdateElemInfo(problem)
  enddo

  call DumpUltra(problem) ! One last dump

  call datastore_delete(ds)

contains

!*************************************************************************
! Subroutine:  GetUserInput
! Purpose   :  Ask for control and output information
!*************************************************************************

subroutine GetUserInput(problem)

  type(datagroup), intent(IN) :: problem
  real(8) :: pratio, dratio

  !********************************!
  ! Get mesh info, and create mesh !
  !********************************!
  print *, "How many zones for the 1D shock tube? "
  !--//read *, numElems
  numElems = 10

  numElems = numElems + 2 ! add an inflow and outflow zone
  numFaces = numElems - 1

!    *(problem->CreateDataObject("numElems")->SetType<int>()->SetDataShape(1)->Allocate()->GetData<int*>()) = numElems
!    *(problem->CreateDataObject("numFaces")->SetType<int>()->SetDataShape(1)->Allocate()->GetData<int*>()) = numFaces

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


!    *(problem->CreateDataObject("pressureRatio")->SetType<double>()->SetDataShape(1)->Allocate()->GetData<double*>()) = pratio
!    *(problem->CreateDataObject("densityRatio")->SetType<double>()->SetDataShape(1)->Allocate()->GetData<double*>()) = dratio
  pressureRatio = pratio
  densityRatio = dratio

  !******************!
  ! Get output  info !
  !******************!
  print *, "How many Ultra dumps would you like? "
!--//    read *, numUltraDumps
  numUltraDumps = 10

  print *, "How many cycles per Ultra dump would you like? "
!--//    read *, numCyclesPerDump
  numCyclesPerDump = 10
!    *(problem->CreateDataObject("numUltraDumps")->SetType<int>()->SetDataShape(1)->Allocate()->GetData<int*>()) = numUltraDumps
!    *(problem->CreateDataObject("numCyclesPerDump")->SetType<int>()->SetDataShape(1)->Allocate()->GetData<int*>()) = numCyclesPerDump
!    *(problem->CreateDataObject("numTotalCycles")->SetType<int>()->SetDataShape(1)->Allocate()->GetData<int*>()) = numUltraDumps * numCyclesPerDump
  numTotalCycles = numUltraDumps * numCyclesPerDump

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

subroutine CreateShockTubeMesh(problem)
  type(datagroup), intent(IN) :: problem
  integer i, k
!  integer numElems = *(prob->GetDataObject("numElems")->GetData<int*>())
!  integer numFaces = *(prob->GetDataObject("numFaces")->GetData<int*>())
  integer inflow(1)
  integer outflow(1)

  ! create element and face classes

!  DataGroup* const elem = prob->CreateDataGroup("elem")->SetDataShape( DataStoreNS::DataShape(numElems))
!  DataGroup* const face = prob->CreateDataGroup("face")->SetDataShape( DataStoreNS::DataShape(numFaces))

  ! set up some important views

  inflow(1) = 1 ! identify inflow elements
!--//  elem->viewCreate("inflow", new IndexSet(1, inflow))
!  elem->CreateDataGroup("inflow")

  outflow(1) = numElems ! identify outflow elements
!--//  elem->viewCreate("outflow", new IndexSet(1, outflow))
!  elem->CreateDataGroup("outflow")

  ! identify shock tube elements - set up basic map
!--//  View *tube = elem->viewCreate("tube", new IndexSet(numElems - 2))
  ! Shift IndexSet indices to identify shocktube elements
  ! (shocktube element numbers are one through numElems-1 inclusive)
!--//  tube->indexSet()->shift(1)
  numTubeElems = numElems - 2
!  DataGroup* const tube = elem->CreateDataGroup("tube")->SetDataShape(DataStoreNS::DataShape(numTubeElems))
!  int* const mapToElems = tube->CreateDataObject("mapToElems")
!                              ->SetType<int>()
!                              ->Allocate()->GetData<int*>()
  allocate(mapToElems(numTubeElems))

  do k=1, numTubeElems
    mapToElems(k) = k + 1
  enddo

  ! Set up some important data relations

  ! Each face connects to two elements
!--//  Relation &faceToElem = *face->relationCreate("faceToElem", 2)
!  std::size_t dims(2) = { numFaces, 2 }
!  DataStoreNS::DataShape desc(2, dims)

!  integer * const faceToElem = face->CreateDataObject("faceToElem")
!                                ->SetDataShape(desc)
!                                ->SetType<int>()
!                                ->Allocate()->GetData<int*>()
  allocate(faceToElem(2, numFaces))

  do i=1, numFaces
    faceToElem(IUPWIND, i) = i
    faceToElem(IDOWNWIND, i) = i + 1
  enddo

  ! Each element connects to two faces
!--//  Relation &elemToFace = *tube->relationCreate("elemToFace", 2)
!  dims(0) = numElems
!  DataStoreNS::DataShape desc2(2, dims)
!  integer * const elemToFace = tube->CreateDataObject("elemToFace")->SetType<int>()->SetDataShape(desc2)->Allocate()->GetData<int*>()
  allocate(elemToFace(2, numElems))

  do i=1, numElems
    elemToFace(IUPWIND, i) = i ! same map as above by coincidence
    elemToFace(IDOWNWIND, i) = i + 1
  enddo

  return
end subroutine CreateShockTubeMesh

!************************************************************************
! Subroutine:  InitializeShockTube
! Purpose   :  Populate the mesh with values
!************************************************************************

subroutine InitializeShockTube(problem)
  type(datagroup), intent(IN) :: problem
  integer i

  integer startTube
  integer endTube
  integer midTube
  real(8) massInitial
  real(8) momentumInitial
  real(8) pressureInitial
  real(8) energyInitial
  real(8) dratio, pratio

  ! These were created in GetUserInput()
!  DataGroup* const elem = (prob->GetDataGroup("elem"))
!  DataGroup* const face = (prob->GetDataGroup("face"))

  ! Create element centered quantities
!  double *mass = elem->CreateDataObject("mass")->SetType<double>()->Allocate()->GetData<double*>()
!  double *momentum = elem->CreateDataObject("momentum")->SetType<double>()->Allocate()->GetData<double*>()
!  double *energy = elem->CreateDataObject("energy")->SetType<double>()->Allocate()->GetData<double*>()
!  double *pressure = elem->CreateDataObject("pressure")->SetType<double>()->Allocate()->GetData<double*>()
  allocate(mass(numElems))
  allocate(momentum(numElems))
  allocate(energy(numElems))
  allocate(pressure(numElems))

  ! Create face centered quantities
!  face->CreateDataObject("F0")->SetType<double>()->Allocate()
!  face->CreateDataObject("F1")->SetType<double>()->Allocate()
!  face->CreateDataObject("F2")->SetType<double>()->Allocate()
  allocate(F0(numFaces))
  allocate(F1(numFaces))
  allocate(F2(numFaces))

!--//  face->fieldCreateReal("F", 3) ! mv, mv^2+P, and v(E+P)

  ! Fill left half with high pressure, right half with low pressure
  startTube = 1
  endTube = numElems  ! elem->GetDataShape().m_dimensions(0)
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
!  dratio = *(prob->GetDataObject("densityRatio")->GetData<real(8)*>())
!  pratio = *(prob->GetDataObject("pressureRatio")->GetData<real(8)*>())
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
!  *(prob->CreateDataObject("time")->SetType<real(8)>()->SetDataShape(1)->Allocate()->GetData<real(8)*>()) = 0.0
!  *(prob->CreateDataObject("cycle")->SetType<int>()->SetDataShape(1)->Allocate()->GetData<int*>()) = 0
!  *(prob->CreateDataObject("dx")->SetType<real(8)>()->SetDataShape(1)->Allocate()->GetData<real(8)*>()) = (1.0 / ((real(8)) endTube))
!  real(8) dx = *(prob->GetDataObject("dx")->GetData<real(8)*>())
!  *(prob->CreateDataObject("dt")->SetType<real(8)>()->SetDataShape(1)->Allocate()->GetData<real(8)*>()) = 0.4 * dx

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

subroutine ComputeFaceInfo(problem)
  type(datagroup), intent(IN) :: problem

  integer i
  integer upWind, downWind, contributor
  real(8) massf, momentumf, energyf, pressuref, c, v, ev, cLocal
!  DataGroup* const face = problem->GetDataGroup("face")
!--//  Relation &faceToElem = *face->relation("faceToElem")
!  integer const * const faceToElem = face->GetDataObject("faceToElem")->GetData<int*>()

!  real(8) * const F0 = face->GetDataObject("F0")->GetData<real(8)*>()
!  real(8) * const F1 = face->GetDataObject("F1")->GetData<real(8)*>()
!  real(8) * const F2 = face->GetDataObject("F2")->GetData<real(8)*>()
!  integer numFaces = face->GetDataShape().m_dimensions(0)

!  DataGroup* const elem = problem->GetDataGroup("elem")
!  real(8) *mass = elem->GetDataObject("mass")->GetData<real(8)*>()
!  real(8) *momentum = elem->GetDataObject("momentum")->GetData<real(8)*>()
!  real(8) *energy = elem->GetDataObject("energy")->GetData<real(8)*>()


  do i=1, numFaces
    ! each face has an upwind and downwind element.
    upWind = faceToElem(IUPWIND, i) ! upwind element
    downWind = faceToElem(IDOWNWIND, i) ! downwind element

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

subroutine UpdateElemInfo(problem)
  type(datagroup), intent(IN) :: problem

  integer i
  integer elemIdx, upWind, downWind

  ! get the element quantities we want to update
!  DataGroup* const elem = problem->GetDataGroup("elem")
!  real(8) * const mass = elem->GetDataObject("mass")->GetData<real(8)*>()
!  real(8) * const momentum = elem->GetDataObject("momentum")->GetData<real(8)*>()
!  real(8) * const energy = elem->GetDataObject("energy")->GetData<real(8)*>()
!  real(8) * const pressure = elem->GetDataObject("pressure")->GetData<real(8)*>()

  ! focus on just the elements within the shock tube
!  DataGroup* const tube = elem->GetDataGroup("tube")
!  integer * const elemToFace = tube->GetDataObject("elemToFace")->GetData<int*>()

!--//  Relation &elemToFace = *tube->relation("elemToFace")
!  integer numTubeElems = tube->GetDataShape().m_dimensions(0)

!--//  int *is = tube->map()
!  integer* const is = tube->GetDataObject("mapToElems")->GetData<int*>()



  ! The element update is calculated as the flux between faces
!  DataGroup* const face = problem->GetDataGroup("face")
!  real(8) *F0 = face->GetDataObject("F0")->GetData<real(8)*>()
!  real(8) *F1 = face->GetDataObject("F1")->GetData<real(8)*>()
!  real(8) *F2 = face->GetDataObject("F2")->GetData<real(8)*>()

!  real(8) dx = *(problem->GetDataObject("dx")->GetData<real(8)*>())
!  real(8) dt = *(problem->GetDataObject("dt")->GetData<real(8)*>())
!  real(8)& time = *(problem->GetDataObject("time")->GetData<real(8)*>())



  do i=1, numTubeElems
     ! recalculate elements in the shocktube, don't touch inflow/outflow
     elemIdx = mapToElems(i)
     ! each element inside the tube has an upwind and downwind face
     upWind = elemToFace(IUPWIND, i) ! upwind face
     downWind = elemToFace(IDOWNWIND, i) ! downwind face

     mass(elemIdx) = mass(elemIdx) - (gammaaInverse * (F0(downWind) - F0(upWind)) * dt / dx)
     momentum(elemIdx) = momentum(elemIdx) - (gammaaInverse * (F1(downWind) - F1(upWind)) * dt / dx)
     energy(elemIdx) = energy(elemIdx) - (gammaaInverse * (F2(downWind) - F2(upWind)) * dt / dx)
     pressure(elemIdx) = (gammaa - 1.0d0) * (energy(elemIdx) - 0.5d0 * momentum(elemIdx) * momentum(elemIdx) / mass(elemIdx))
  enddo

  ! update the time
  time = time + dt

  return
end subroutine UpdateElemInfo

subroutine DumpUltra( problem )
  type(datagroup), intent(IN) :: problem
   integer fp
   character(100) fname

!--#if 0
!--   char *tail
!--
!--!--//   VHashTraverse_t content
!--
!--!   const DataGroup* const elem = prob->GetDataGroup("elem")
!--
!--   fname = "problem"
!--
!--   ! Skip past the junk
!--!   for (tail=fname isalpha(*tail) ++tail)
!--
!--!   sprintf(tail, "_%04d", *(prob->GetDataObject("cycle")->GetData<int*>()) )
!--
!--   
!--   fp = open(fname, "w")) == NULL)
!--   {
!--      print *, "Could not open file " // trim(fname) // ". Aborting.")
!--      stop
!--   }
!--
!--   write(fp, "# Ultra Plot File\n")
!--   write(fp, "# Problem: %s\n", "problem" )
!--
!--
!--
!--   {
!--   const DataGroup::dataArrayType& dataObjects = prob->GetDataObjects()
!--   const DataGroup::lookupType& dataObjectLookup = prob->GetDataObjectLookup()
!--
!--   DataGroup::dataArrayType::const_iterator obj=dataObjects.begin()
!--   DataGroup::lookupType::const_iterator lookup=dataObjectLookup.begin()
!--
!--   for( obj!=dataObjects.end() ++obj, ++lookup )
!--   {
!--     const int length = (*obj)->GetDataShape().m_dimensions(0)
!--     const std::string& name = lookup->first
!--     if( length <= 1 )
!--     {
!--       if( (*obj)->GetType() == DataStoreNS::rtTypes::int32_id )
!--       {
!--         write(fp, "# %s = %d\n", name.c_str(), *((*obj)->GetData<int*>()))
!--       }
!--       else if( (*obj)->GetType() == DataStoreNS::rtTypes::real64_id )
!--       {
!--         write(fp, "# %s = %f\n", name.c_str(), *((*obj)->GetData<real(8)*>()))
!--       }
!--     }
!--   }
!--   }
!--
!--   {
!--!--//   for( auto obj : elem->GetDataObjects() )
!--   const DataGroup::dataArrayType& dataObjects = elem->GetDataObjects()
!--   const DataGroup::lookupType& dataObjectLookup = elem->GetDataObjectLookup()
!--
!--   DataGroup::dataArrayType::const_iterator obj=dataObjects.begin()
!--   DataGroup::lookupType::const_iterator lookup=dataObjectLookup.begin()
!--
!--   for( obj!=dataObjects.end() ++obj, ++lookup )
!--   {
!--     const int length = (*obj)->GetDataShape().m_dimensions(0)
!--     const std::string& name = lookup->first
!--     fprintf(fp, "# %s\n", name.c_str() )
!--     if( (*obj)->GetType() == DataStoreNS::rtTypes::int32_id )
!--     {
!--       int const * const data = (*obj)->GetData<int*>()
!--       for ( int i=0 i<length ++i)
!--          write(fp, "%f %f\n", (real(8)) i, (real(8)) data(i))
!--       write(fp, "\n")
!--     }
!--     else if( (*obj)->GetType() == DataStoreNS::rtTypes::real64_id )
!--     {
!--       real(8) const * const data = (*obj)->GetData<real(8)*>()
!--       for ( int i=0 i<length ++i)
!--          write(fp, "%f %f\n", (real(8)) i, (real(8)) data(i) )
!--       write(fp, "\n")
!--     }
!--   }
!--   }
!--
!--   close(fp)
!--#endif
   return
 end subroutine DumpUltra

end program main
