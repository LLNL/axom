/**************************************************************************
 * Program:  ShockTube.C
 * Purpose:  1D shock tube, split flux Euler equations
 *
 *         | m  |            |    mv    |
 *     Q = | mv |        F = | mv^2 + P |
 *         | E  |            |  v(E+P)  |
 *
 *     P = (gamma - 1.0)[E - 0.5 mv^2 ]
 *
 *             Cp
 *     gamma = --     m = mass/volume   v = velocity
 *             Cv
 *
 *     All quantiites are non-dimensionalized.
 *
 *     @Q   @F    @Q   @F @Q
 *     -- + -- =  -- + -- -- = 0
 *     @t   @x    @t   @Q @x
 *
 *************************************************************************/

//#include "Vista.h"
//#include "View.h"
#include "../src/Types.hpp"
#include "../src/DataStore.hpp"

#include <stdio.h>
#include <math.h>
#include "stdlib.h"
#define UPWIND   0
#define DOWNWIND 1

const double gammaa = M_SQRT2;
const double gammaaInverse = M_SQRT1_2;

using DataStoreNS::DataView;
using DataStoreNS::DataGroup;

/**************************************************************************
 * Subroutine:  GetUserInput
 * Purpose   :  Ask for control and output information
 *************************************************************************/

void GetUserInput(DataStoreNS::DataGroup* const problem)
{
  /**********************************/
  /* Get mesh info, and create mesh */
  /**********************************/
  {
    int numElems, numFaces;

    printf("How many zones for the 1D shock tube? ");
    //scanf("%d", &numElems);
    numElems = 10;

    numElems += 2; /* add an inflow and outflow zone */
    numFaces = numElems - 1;

    *(problem->CreateDataView("numElems")->SetType<int>()->SetDataShape(1)->Allocate()->GetData<int*>()) = numElems;
    *(problem->CreateDataView("numFaces")->SetType<int>()->SetDataShape(1)->Allocate()->GetData<int*>()) = numFaces;
  }

  /********************/
  /* Get pyhsics info */
  /********************/
  {
    double pratio = -1.0, dratio = -1.0;

#if 0
    printf("What cfl number would you like to use? ");
    scanf("%lf", &cfl);
#endif

    while (pratio < 0.0 || pratio > 1.0)
    {
      printf("What pressure ratio would you like (0 <= x <= 1)? ");
      //scanf("%lf", &pratio);
      pratio = 0.5;
    }

    while (dratio < 0.0 || dratio > 1.0)
    {
      printf("What density ratio would you like (0 <= x <= 1)? ");
      //scanf("%lf", &dratio);
      dratio = 0.5;
    }


    *(problem->CreateDataView("pressureRatio")->SetType<double>()->SetDataShape(1)->Allocate()->GetData<double*>()) = pratio;
    *(problem->CreateDataView("densityRatio")->SetType<double>()->SetDataShape(1)->Allocate()->GetData<double*>()) = dratio;


  }

  /********************/
  /* Get output  info */
  /********************/
  {
    int numUltraDumps, numCyclesPerDump;

    printf("How many Ultra dumps would you like? ");
//    scanf("%d", &numUltraDumps);
    numUltraDumps = 10;

    printf("How many cycles per Ultra dump would you like? ");
//    scanf("%d", &numCyclesPerDump);
    numCyclesPerDump = 10;
    *(problem->CreateDataView("numUltraDumps")->SetType<int>()->SetDataShape(1)->Allocate()->GetData<int*>()) = numUltraDumps;
    *(problem->CreateDataView("numCyclesPerDump")->SetType<int>()->SetDataShape(1)->Allocate()->GetData<int*>()) = numCyclesPerDump;
    *(problem->CreateDataView("numTotalCycles")->SetType<int>()->SetDataShape(1)->Allocate()->GetData<int*>()) = numUltraDumps * numCyclesPerDump;

  }

  return;
}

/**************************************************************************
 * Subroutine:  CreateShockTubeMesh
 * Purpose   :  Build an empty mesh for the shock tube


 Gaps between elements are faces
 |
 -------------------------------
 |   |   |             |   |   |

 ### ### ### ###       ### ### ### ###
 ### ### ### ###  ...  ### ### ### ###  <--- 1D Shock tube model
 ### ### ### ###       ### ### ### ###

 |  |                           |  |
 |  -----------------------------  |
 Inflow           |               Outflow
 Element      Tube Elements       Element


 *************************************************************************/

void CreateShockTubeMesh(DataGroup * const prob)
{
  int i;
  int numElems = *(prob->GetDataView("numElems")->GetData<int*>());
  int numFaces = *(prob->GetDataView("numFaces")->GetData<int*>());
  int inflow[1];
  int outflow[1];

  /* create element and face classes */

  DataGroup* const elem = prob->CreateDataGroup("elem")->SetDataShape( DataStoreNS::DataShape(numElems));
  DataGroup* const face = prob->CreateDataGroup("face")->SetDataShape( DataStoreNS::DataShape(numFaces));;

  /* set up some important views */

  inflow[0] = 0; /* identify inflow elements */
//  elem->viewCreate("inflow", new IndexSet(1, inflow));
  elem->CreateDataGroup("inflow");

  outflow[0] = numElems - 1; /* identify outflow elements */
//  elem->viewCreate("outflow", new IndexSet(1, outflow));
  elem->CreateDataGroup("outflow");

  /* identify shock tube elements - set up basic map */
//  View *tube = elem->viewCreate("tube", new IndexSet(numElems - 2));
  /* Shift IndexSet indices to identify shocktube elements */
  /* (shocktube element numbers are one through numElems-1 inclusive) */
//  tube->indexSet()->shift(1);
  std::size_t numTubeElems = (numElems - 2);
  DataGroup* const tube = elem->CreateDataGroup("tube")->SetDataShape(DataStoreNS::DataShape(numTubeElems));
  int* const mapToElems = tube->CreateDataView("mapToElems")
                              ->SetType<int>()
                              ->Allocate()->GetData<int*>();
  for ( unsigned int k = 0u; k < numTubeElems; ++k)
  {
    mapToElems[k] = k + 1;
  }

  /* Set up some important data relations */

  /* Each face connects to two elements */
//  Relation &faceToElem = *face->relationCreate("faceToElem", 2);
  std::size_t dims[2] = { numFaces, 2 };
  DataStoreNS::DataShape desc(2, dims);

  int * const faceToElem = face->CreateDataView("faceToElem")
                                ->SetDataShape(desc)
                                ->SetType<int>()
                                ->Allocate()->GetData<int*>();

  for (i = 0; i < numFaces; ++i)
  {
    faceToElem[i * 2 + UPWIND] = i;
    faceToElem[i * 2 + DOWNWIND] = i + 1;
  }

  /* Each element connects to two faces */ //
//  Relation &elemToFace = *tube->relationCreate("elemToFace", 2);
  dims[0] = numElems;
  DataStoreNS::DataShape desc2(2, dims);
  int * const elemToFace = tube->CreateDataView("elemToFace")->SetType<int>()->SetDataShape(desc2)->Allocate()->GetData<int*>();

  for (i = 0; i < numElems; ++i)
  {
    elemToFace[i * 2 + UPWIND] = i; /* same map as above by coincidence */
    elemToFace[i * 2 + DOWNWIND] = i + 1;
  }

  return;
}

/**************************************************************************
 * Subroutine:  InitializeShockTube
 * Purpose   :  Populate the mesh with values
 *************************************************************************/

void InitializeShockTube(DataGroup * const prob)
{
  int i;

  /* These were created in GetUserInput() */
  DataGroup* const elem = (prob->GetDataGroup("elem"));
  DataGroup* const face = (prob->GetDataGroup("face"));

  /* Create element centered quantities */
  double *mass = elem->CreateDataView("mass")->SetType<double>()->Allocate()->GetData<double*>();
  double *momentum = elem->CreateDataView("momentum")->SetType<double>()->Allocate()->GetData<double*>();
  double *energy = elem->CreateDataView("energy")->SetType<double>()->Allocate()->GetData<double*>();
  double *pressure = elem->CreateDataView("pressure")->SetType<double>()->Allocate()->GetData<double*>();

  /* Create face centered quantities */
  face->CreateDataView("F0")->SetType<double>()->Allocate();
  face->CreateDataView("F1")->SetType<double>()->Allocate();
  face->CreateDataView("F2")->SetType<double>()->Allocate();

//  face->fieldCreateReal("F", 3); /* mv, mv^2+P, and v(E+P) */

  /* Fill left half with high pressure, right half with low pressure */
  int startTube = 0;
  int endTube = elem->GetDataShape().m_dimensions[0];
  int midTube = endTube / 2;

  /* Non-dimensionalized reference values */
  double massInitial = 1.0;
  double momentumInitial = 0.0;
  double pressureInitial = gammaaInverse;
  double energyInitial = pressureInitial / (gammaa - 1.0);

  /* Initialize zonal quantities*/
  for (i = startTube; i < midTube; ++i)
  {
    mass[i] = massInitial;
    momentum[i] = momentumInitial;
    pressure[i] = pressureInitial;
    energy[i] = energyInitial;
  }

  /* adjust parameters for low pressure portion of tube */
  double dratio = *(prob->GetDataView("densityRatio")->GetData<double*>());
  double pratio = *(prob->GetDataView("pressureRatio")->GetData<double*>());

  massInitial *= dratio;
  pressureInitial *= pratio;
  energyInitial = pressureInitial / (gammaa - 1.0);

  for (i = midTube; i < endTube; ++i)
  {
    mass[i] = massInitial;
    momentum[i] = momentumInitial;
    pressure[i] = pressureInitial;
    energy[i] = energyInitial;
  }

  /* Create needed time info */
  *(prob->CreateDataView("time")->SetType<double>()->SetDataShape(1)->Allocate()->GetData<double*>()) = 0.0;
  *(prob->CreateDataView("cycle")->SetType<int>()->SetDataShape(1)->Allocate()->GetData<int*>()) = 0;
  *(prob->CreateDataView("dx")->SetType<double>()->SetDataShape(1)->Allocate()->GetData<double*>()) = (1.0 / ((double) endTube));
  double dx = *(prob->GetDataView("dx")->GetData<double*>());
  *(prob->CreateDataView("dt")->SetType<double>()->SetDataShape(1)->Allocate()->GetData<double*>()) = 0.4 * dx;



  return;
}

/**************************************************************************
 * Subroutine:  ComputeFaceInfo
 * Purpose   :  Compute F quantities at faces.
 *
 *  @F   @F0   @F1   @F2
 *  -- = --- + --- + ---
 *  @x   @x    @x    @x
 *
 *  Calculate F0, F1 and F2 at the face centers.
 *
 *************************************************************************/

void ComputeFaceInfo(DataGroup * const problem)
{
  int i;
  DataGroup* const face = problem->GetDataGroup("face");
//  Relation &faceToElem = *face->relation("faceToElem");
  int const * const faceToElem = face->GetDataView("faceToElem")->GetData<int*>();

  double * const F0 = face->GetDataView("F0")->GetData<double*>();
  double * const F1 = face->GetDataView("F1")->GetData<double*>();
  double * const F2 = face->GetDataView("F2")->GetData<double*>();
  int numFaces = face->GetDataShape().m_dimensions[0];

  DataGroup* const elem = problem->GetDataGroup("elem");
  double *mass = elem->GetDataView("mass")->GetData<double*>();
  double *momentum = elem->GetDataView("momentum")->GetData<double*>();
  double *energy = elem->GetDataView("energy")->GetData<double*>();

  for (i = 0; i < numFaces; ++i)
  {
    /* each face has an upwind and downwind element. */
    int upWind = faceToElem[i * 2 + UPWIND]; /* upwind element */
    int downWind = faceToElem[i * 2 + DOWNWIND]; /* downwind element */

    /* calculate face centered quantities */
    double massf = 0.5 * (mass[upWind] + mass[downWind]);
    double momentumf = 0.5 * (momentum[upWind] + momentum[downWind]);
    double energyf = 0.5 * (energy[upWind] + energy[downWind]);
    double pressuref = (gammaa - 1.0) * (energyf - 0.5 * momentumf * momentumf / massf);
    double c = sqrt(gammaa * pressuref / massf);
    double v = momentumf / massf;
    double ev;
    double cLocal;
    int contributor;

    /* Now that we have the wave speeds, we might want to */
    /* look for the max wave speed here, and update dt */
    /* appropriately right before leaving this function. */
    /* ... */

    /* OK, calculate face quantities */

    F0[i] = F1[i] = F2[i] = 0.0;

    contributor = ((v >= 0.0) ? upWind : downWind);
    massf = mass[contributor];
    momentumf = momentum[contributor];
    energyf = energy[contributor];
    pressuref = energyf - 0.5 * momentumf * momentumf / massf;
    ev = v * (gammaa - 1.0);

    F0[i] += ev * massf;
    F1[i] += ev * momentumf;
    F2[i] += ev * (energyf - pressuref);

    contributor = ((v + c >= 0.0) ? upWind : downWind);
    massf = mass[contributor];
    momentumf = momentum[contributor];
    energyf = energy[contributor];
    pressuref = (gammaa - 1.0) * (energyf - 0.5 * momentumf * momentumf / massf);
    ev = 0.5 * (v + c);
    cLocal = sqrt(gammaa * pressuref / massf);

    F0[i] += ev * massf;
    F1[i] += ev * (momentumf + massf * cLocal);
    F2[i] += ev * (energyf + pressuref + momentumf * cLocal);

    contributor = ((v - c >= 0.0) ? upWind : downWind);
    massf = mass[contributor];
    momentumf = momentum[contributor];
    energyf = energy[contributor];
    pressuref = (gammaa - 1.0) * (energyf - 0.5 * momentumf * momentumf / massf);
    ev = 0.5 * (v - c);
    cLocal = sqrt(gammaa * pressuref / massf);

    F0[i] += ev * massf;
    F1[i] += ev * (momentumf - massf * cLocal);
    F2[i] += ev * (energyf + pressuref - momentumf * cLocal);
  }

  return;
}

/**************************************************************************
 * Subroutine:  UpdateElemInfo
 * Purpose   :  Q(elem) = Q(elem) + deltaQ(elem)
 *
 *  deltaQ(elem) = - (F(downWindFace) - F(upWindFace)) * dt / dx ;
 *
 *************************************************************************/

void UpdateElemInfo(DataGroup * const problem)
{
  int i;

  /* get the element quantities we want to update */
  DataGroup* const elem = problem->GetDataGroup("elem");
  double * const mass = elem->GetDataView("mass")->GetData<double*>();
  double * const momentum = elem->GetDataView("momentum")->GetData<double*>();
  double * const energy = elem->GetDataView("energy")->GetData<double*>();
  double * const pressure = elem->GetDataView("pressure")->GetData<double*>();

  /* focus on just the elements within the shock tube */
  DataGroup* const tube = elem->GetDataGroup("tube");
  int * const elemToFace = tube->GetDataView("elemToFace")->GetData<int*>();

//  Relation &elemToFace = *tube->relation("elemToFace");
  int numTubeElems = tube->GetDataShape().m_dimensions[0];

//  int *is = tube->map();
  int* const is = tube->GetDataView("mapToElems")->GetData<int*>();



  /* The element update is calculated as the flux between faces */
  DataGroup* const face = problem->GetDataGroup("face");
  double *F0 = face->GetDataView("F0")->GetData<double*>();
  double *F1 = face->GetDataView("F1")->GetData<double*>();
  double *F2 = face->GetDataView("F2")->GetData<double*>();

  double dx = *(problem->GetDataView("dx")->GetData<double*>());
  double dt = *(problem->GetDataView("dt")->GetData<double*>());
  double& time = *(problem->GetDataView("time")->GetData<double*>());

  for (i = 0; i < numTubeElems; ++i)
  {
    /* recalculate elements in the shocktube, don't touch inflow/outflow */
    int elemIdx = is[i];

    /* each element inside the tube has an upwind and downwind face */
    int upWind = elemToFace[i * 2 + UPWIND]; /* upwind face */
    int downWind = elemToFace[i * 2 + DOWNWIND]; /* downwind face */

    mass[elemIdx] -= gammaaInverse * (F0[downWind] - F0[upWind]) * dt / dx;
    momentum[elemIdx] -= gammaaInverse * (F1[downWind] - F1[upWind]) * dt / dx;
    energy[elemIdx] -= gammaaInverse * (F2[downWind] - F2[upWind]) * dt / dx;
    pressure[elemIdx] = (gammaa - 1.0) * (energy[elemIdx] - 0.5 * momentum[elemIdx] * momentum[elemIdx] / mass[elemIdx]);
  }

  /* update the time */
  time += dt;

  return;
}

#include <stdio.h>
#include <string.h>
#include <ctype.h>


void DumpUltra( const DataGroup * const prob)
{
   FILE *fp ;
   char fname[100] ;
   char *tail ;

//   VHashTraverse_t content ;

   const DataGroup* const elem = prob->GetDataGroup("elem") ;

   strcpy(fname, "problem" ) ;

   /* Skip past the junk */
   for (tail=fname; isalpha(*tail); ++tail) ;

   sprintf(tail, "_%04d", *(prob->GetDataView("cycle")->GetData<int*>()) ) ;

   if ((fp = fopen(fname, "w")) == NULL)
   {
      printf("Could not open file %s. Aborting.\n", fname) ;
      exit (-1) ;
   }

   fprintf(fp, "# Ultra Plot File\n") ;
   fprintf(fp, "# Problem: %s\n", "problem" ) ;



   {
   const DataGroup::dataArrayType& dataObjects = prob->GetDataViews();
   const DataGroup::lookupType& dataObjectLookup = prob->GetDataViewLookup();

   DataGroup::dataArrayType::const_iterator obj=dataObjects.begin();
   DataGroup::lookupType::const_iterator lookup=dataObjectLookup.begin();

   for(  ; obj!=dataObjects.end() ; ++obj, ++lookup )
   {
     const int length = (*obj)->GetDataShape().m_dimensions[0];
     const std::string& name = lookup->first;
     if( length <= 1 )
     {
       if( (*obj)->GetType() == DataStoreNS::rtTypes::int32_id )
       {
         fprintf(fp, "# %s = %d\n", name.c_str(), *((*obj)->GetData<int*>())) ;
       }
       else if( (*obj)->GetType() == DataStoreNS::rtTypes::real64_id )
       {
         fprintf(fp, "# %s = %f\n", name.c_str(), *((*obj)->GetData<double*>())) ;
       }
     }
   }
   }

   {
//   for( auto obj : elem->GetDataViews() )
   const DataGroup::dataArrayType& dataObjects = elem->GetDataViews();
   const DataGroup::lookupType& dataObjectLookup = elem->GetDataViewLookup();

   DataGroup::dataArrayType::const_iterator obj=dataObjects.begin();
   DataGroup::lookupType::const_iterator lookup=dataObjectLookup.begin();

   for(  ; obj!=dataObjects.end() ; ++obj, ++lookup )
   {
     const int length = (*obj)->GetDataShape().m_dimensions[0];
     const std::string& name = lookup->first;
     fprintf(fp, "# %s\n", name.c_str() ) ;
     if( (*obj)->GetType() == DataStoreNS::rtTypes::int32_id )
     {
       int const * const data = (*obj)->GetData<int*>();
       for ( int i=0; i<length; ++i)
          fprintf(fp, "%f %f\n", (double) i, (double) data[i]) ;
       fprintf(fp, "\n") ;
     }
     else if( (*obj)->GetType() == DataStoreNS::rtTypes::real64_id )
     {
       double const * const data = (*obj)->GetData<double*>();
       for ( int i=0; i<length; ++i)
          fprintf(fp, "%f %f\n", (double) i, (double) data[i]) ;
       fprintf(fp, "\n") ;
     }
   }
   }

   fclose(fp) ;

   return ;
}
/**************************************************************************
 * Subroutine:  main
 * Purpose   :  Simulate a 1D Shock Tube using split flux Euler formulation
 *************************************************************************/

int main(void)
{
 // extern void DumpUltra(DataGroup * const prob);

  /* We should be able to parallelize pretty easily by */
  /* adding an MPI_Init() here, modifying the setup slightly, */
  /* adding a communication subroutine in the main loop, */
  /* and calling MPI_Finalize() at the end of main() */

  DataStoreNS::DataStore DATASTORE;
  DataStoreNS::DataStore* const dataStore = &DATASTORE;
  DataStoreNS::DataGroup* const rootGroup = dataStore->GetRootDataGroup();
  DataStoreNS::DataGroup* const problem = rootGroup->CreateDataGroup("problem");

  GetUserInput(problem);
  CreateShockTubeMesh(problem);
  InitializeShockTube(problem);


  /* use a reference when you want to update the param directly */
  int* const currCycle = problem->GetDataView("cycle")->GetData<int*>();

  int numTotalCycles = *(problem->GetDataView("numTotalCycles")->GetData<int*>());
  int dumpInterval = *(problem->GetDataView("numCyclesPerDump")->GetData<int*>());

  for (*currCycle = 0; *currCycle < numTotalCycles; ++(*currCycle) )
  {
    /* dump the ultra file, based on the user chosen attribute mask */
    if ( (*currCycle) % dumpInterval == 0)
      DumpUltra(problem);

    ComputeFaceInfo(problem);
    UpdateElemInfo(problem);
  }

  DumpUltra(problem); /* One last dump */

  return 0;
}

