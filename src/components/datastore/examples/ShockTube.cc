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
#include "../src/DataBuffer.hpp"

#include <stdio.h>
#include <math.h>
#include "stdlib.h"
#define UPWIND   0
#define DOWNWIND 1

const double gammaa = M_SQRT2;
const double gammaaInverse = M_SQRT1_2;

using sidre::DataView;
using sidre::DataGroup;
using sidre::DataBuffer;
using namespace conduit;


void CreateScalarIntBufferViewAndSetVal( DataGroup* const grp, const std::string& name, int32 const value )
{
  DataBuffer* const buffer = grp->GetDataStore()->CreateBuffer()
                                                ->Declare(DataType::int32())
                                                ->Allocate();


    DataView* const view = grp->CreateView(name, buffer)->Apply(DataType::int32());
  *(view->GetNode().as_int32_ptr()) = value;

}


void CreateScalarFloatBufferViewAndSetVal( DataGroup* const grp, const std::string& name, float64 const value )
{
  DataBuffer* const buffer = grp->GetDataStore()->CreateBuffer()
                                                ->Declare(DataType::float64())
                                                ->Allocate();


  DataView* const view = grp->CreateView(name,buffer)->Apply(DataType::float64());
  *(view->GetNode().as_float64_ptr()) = value;

}




/**************************************************************************
 * Subroutine:  GetUserInput
 * Purpose   :  Ask for control and output information
 *************************************************************************/

void GetUserInput(sidre::DataGroup* const problem)
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

    // create buffer and view, and set value
    CreateScalarIntBufferViewAndSetVal( problem, "numElems", numElems );

    CreateScalarIntBufferViewAndSetVal( problem, "numFaces", numFaces );

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


    CreateScalarFloatBufferViewAndSetVal( problem, "pressureRatio", pratio );

    CreateScalarFloatBufferViewAndSetVal( problem, "densityRatio", dratio );


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



    CreateScalarIntBufferViewAndSetVal( problem, "numUltraDumps", numUltraDumps );

    CreateScalarIntBufferViewAndSetVal( problem, "numCyclesPerDump", numCyclesPerDump );

    CreateScalarIntBufferViewAndSetVal( problem, "numTotalCycles", numUltraDumps * numCyclesPerDump );


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
  int32 const numElems = prob->GetView("numElems")->GetNode().as_int32();
  int32 const numFaces = prob->GetView("numFaces")->GetNode().as_int32();

  int inflow[1];
  int outflow[1];

  /* create element and face classes */

  DataGroup* const elem = prob->CreateGroup("elem");
  DataGroup* const face = prob->CreateGroup("face");

  /* set up some important views */

  inflow[0] = 0; /* identify inflow elements */
//  elem->viewCreate("inflow", new IndexSet(1, inflow));
  elem->CreateGroup("inflow");

  outflow[0] = numElems - 1; /* identify outflow elements */
//  elem->viewCreate("outflow", new IndexSet(1, outflow));
  elem->CreateGroup("outflow");



  std::size_t numTubeElems = (numElems - 2);
  DataGroup* const tube = elem->CreateGroup("tube");//->SetDataShape(DataStoreNS::DataShape(numTubeElems));

  DataBuffer* const mapToElemsBuffer = elem->GetDataStore()->CreateBuffer()
                                                           ->Declare(DataType::int32(numTubeElems))
                                                           ->Allocate();


  DataView* const mapToElemsView = tube->CreateView("mapToElems", mapToElemsBuffer)
                                       ->Apply(DataType::int32(numTubeElems));

  uint32* const mapToElems = mapToElemsView->GetNode().as_uint32_ptr();


  for ( unsigned int k = 0u; k < numTubeElems; ++k)
  {
    mapToElems[k] = k + 1;
  }

  /* Set up some important data relations */

  /* Each face connects to two elements */

  DataBuffer* const faceToElemBuffer = face->GetDataStore()->CreateBuffer()
                                                           ->Declare(DataType::int32(2*numFaces,0,8))
                                                           ->Allocate();

  int32* const faceToElem = face->CreateView("faceToElem", faceToElemBuffer)
                                              ->Apply(DataType::int32(2*numFaces,0,8))
                                              ->GetNode().as_int32_ptr();

  for (i = 0; i < numFaces; ++i)
  {
    faceToElem[i * 2 + UPWIND] = i;
    faceToElem[i * 2 + DOWNWIND] = i + 1;
  }

  /* Each element connects to two faces */ //
//  Relation &elemToFace = *tube->relationCreate("elemToFace", 2);
  DataBuffer* const elemToFaceBuffer = tube->GetDataStore()->CreateBuffer()
                                                           ->Declare(DataType::int32(2*numElems,0,8))
                                                           ->Allocate();


  int32* elemToFace = tube->CreateView("elemToFace", faceToElemBuffer)
                                       ->Apply(DataType::int32(2*numElems,0,8))
                                       ->GetNode().as_int32_ptr();

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
  DataGroup* const elem = (prob->GetGroup("elem"));
  DataGroup* const face = (prob->GetGroup("face"));

  /* Create element centered quantities */


  DataBuffer* buffer = nullptr;

  int32 const numElems = prob->GetView("numElems")->GetNode().as_int32();
  int32 const numFaces = prob->GetView("numFaces")->GetNode().as_int32();

  buffer = elem->GetDataStore()->CreateBuffer()
                               ->Declare(DataType::float64(numElems))
                               ->Allocate();


  float64* const mass = elem->CreateView("mass", buffer)
                             ->Apply(DataType::float64(numElems))
                            ->GetNode().as_float64_ptr();


  buffer = elem->GetDataStore()->CreateBuffer()
                               ->Declare(DataType::float64(numElems))
                               ->Allocate();

  float64* const momentum = elem->CreateView("momentum", buffer)
                                ->Apply(DataType::float64(numElems))
                                ->GetNode().as_float64_ptr();


  buffer = elem->GetDataStore()->CreateBuffer()
                               ->Declare(DataType::float64(numElems))
                               ->Allocate();

  float64* const energy = elem->CreateView("energy", buffer)
                              ->Apply(DataType::float64(numElems))
                              ->GetNode().as_float64_ptr();


  buffer = elem->GetDataStore()->CreateBuffer()
                               ->Declare(DataType::float64(numElems))
                               ->Allocate();

  float64* const pressure = elem->CreateView("pressure", buffer)
                                ->Apply(DataType::float64(numElems))
                                ->GetNode().as_float64_ptr();


  /* Create face centered quantities */

  buffer = face->GetDataStore()->CreateBuffer()
                               ->Declare(DataType::float64(numFaces))
                               ->Allocate();

  face->CreateView("F0", buffer)
      ->Apply(DataType::float64(numFaces));

  buffer = face->GetDataStore()->CreateBuffer()
                               ->Declare(DataType::float64(numFaces))
                               ->Allocate();

  face->CreateView("F1", buffer)
       ->Apply(DataType::float64(numFaces));

  buffer = face->GetDataStore()->CreateBuffer()
                               ->Declare(DataType::float64(numFaces))
                               ->Allocate();

  face->CreateView("F2", buffer)
    ->Apply(DataType::float64(numFaces));


//  face->fieldCreateReal("F", 3); /* mv, mv^2+P, and v(E+P) */

  /* Fill left half with high pressure, right half with low pressure */
  int startTube = 0;
  int endTube = numElems;
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
  double dratio = *(prob->GetView("densityRatio")->GetNode().as_float64_ptr());
  double pratio = *(prob->GetView("pressureRatio")->GetNode().as_float64_ptr());

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

  CreateScalarFloatBufferViewAndSetVal( prob, "time", 0.0 );
  CreateScalarIntBufferViewAndSetVal( prob, "cycle", 0 );

  CreateScalarFloatBufferViewAndSetVal( prob, "dx", (1.0 / ((double) endTube)) );
  double dx = *(prob->GetView("dx")->GetNode().as_float64_ptr());
  CreateScalarFloatBufferViewAndSetVal( prob, "dt", 0.4 * dx );



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
  DataGroup* const face = problem->GetGroup("face");
//  Relation &faceToElem = *face->relation("faceToElem");
  int32 const* const faceToElem = face->GetView("faceToElem")->GetNode().as_int32_ptr();

  float64* const F0 = face->GetView("F0")->GetNode().as_float64_ptr();
  float64* const F1 = face->GetView("F1")->GetNode().as_float64_ptr();
  float64* const F2 = face->GetView("F2")->GetNode().as_float64_ptr();
  int numFaces = face->GetView("F0")->GetNode().as_float64_array().number_of_elements();

  DataGroup* const elem = problem->GetGroup("elem");
  float64* const mass = elem->GetView("mass")->GetNode().as_float64_ptr();
  float64* const momentum = elem->GetView("momentum")->GetNode().as_float64_ptr();
  float64* const energy = elem->GetView("energy")->GetNode().as_float64_ptr();

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
  DataGroup* const elem = problem->GetGroup("elem");
  float64 * const mass = elem->GetView("mass")->GetNode().as_float64_ptr();
  float64 * const momentum = elem->GetView("momentum")->GetNode().as_float64_ptr();
  float64 * const energy = elem->GetView("energy")->GetNode().as_float64_ptr();
  float64 * const pressure = elem->GetView("pressure")->GetNode().as_float64_ptr();

  /* focus on just the elements within the shock tube */
  DataGroup* const tube = elem->GetGroup("tube");
  int32 * const elemToFace = tube->GetView("elemToFace")->GetNode().as_int32_ptr();

//  Relation &elemToFace = *tube->relation("elemToFace");
  int numTubeElems = tube->GetView("mapToElems")->GetNode().as_int32_array().number_of_elements();

//  int *is = tube->map();
  int32* const is = tube->GetView("mapToElems")->GetNode().as_int32_ptr();



  /* The element update is calculated as the flux between faces */
  DataGroup* const face = problem->GetGroup("face");
  float64* const F0 = face->GetView("F0")->GetNode().as_float64_ptr();
  float64* const F1 = face->GetView("F1")->GetNode().as_float64_ptr();
  float64* const F2 = face->GetView("F2")->GetNode().as_float64_ptr();

  double const dx = *(problem->GetView("dx")->GetNode().as_float64_ptr());
  double const dt = *(problem->GetView("dt")->GetNode().as_float64_ptr());
  float64& time = *(problem->GetView("time")->GetNode().as_float64_ptr());

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


void DumpUltra( DataGroup * const prob)
{
#if 1
   FILE *fp ;
   char fname[100] ;
   char *tail ;

//   VHashTraverse_t content ;

   DataGroup* const elem = prob->GetGroup("elem") ;

   strcpy(fname, "problem" ) ;

   /* Skip past the junk */
   for (tail=fname; isalpha(*tail); ++tail) ;

   sprintf(tail, "_%04d.ult", *(prob->GetView("cycle")->GetNode().as_int32_ptr()) ) ;

   if ((fp = fopen(fname, "w")) == NULL)
   {
      printf("Could not open file %s. Aborting.\n", fname) ;
      exit (-1) ;
   }

   fprintf(fp, "# Ultra Plot File\n") ;
   fprintf(fp, "# Problem: %s\n", "problem" ) ;


     for(size_t i=0;i<prob->CountViews();i++)
     {
         DataView * const view = prob->GetView(i);
         const int length = view->GetDescriptor().dtype().number_of_elements();
         const std::string& name = view->GetName();
         if( length <= 1 )
         {
             if( view->GetDescriptor().dtype().id() == DataType::INT32_T )
             {
                 fprintf(fp, "# %s = %d\n",
                         name.c_str(),
                         view->GetNode().as_int32());
             }
             else if( view->GetDescriptor().dtype().id() == 
                      DataType::FLOAT64_T )
             {
                 fprintf(fp, "# %s = %f\n",
                         name.c_str(),
                         view->GetNode().as_float64());
            }
        }
    }


    for(size_t i=0;i<elem->CountViews();i++)
    {
         DataView * const view = elem->GetView(i);
         const int length = view->GetDescriptor().dtype().number_of_elements();
         const std::string& name = view->GetName();
         fprintf(fp, "# %s\n", name.c_str() ) ;
         
         if( view->GetDescriptor().dtype().id() == DataType::INT32_T )
         {
             int32 const * const data = view->GetNode().as_int32_ptr();
             for ( int i=0; i<length; ++i)
             {
                 fprintf(fp, "%f %f\n", (double) i, (double) data[i]) ;
             }
             fprintf(fp, "\n") ;
         }
         else if( view->GetDescriptor().dtype().id() == DataType::FLOAT64_T )
         {
             float64 const * const data = view->GetNode().as_float64_ptr();
             for ( int i=0; i<length; ++i)
             {
                 fprintf(fp, "%f %f\n", (double) i, (double) data[i]) ;
             }
             fprintf(fp, "\n") ;
         }
   }


   fclose(fp) ;
#endif
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

  sidre::DataStore DATASTORE;
  sidre::DataStore* const dataStore = &DATASTORE;
  sidre::DataGroup* const rootGroup = dataStore->GetRoot();
  sidre::DataGroup* const problem = rootGroup->CreateGroup("problem");

  GetUserInput(problem);
  CreateShockTubeMesh(problem);
  InitializeShockTube(problem);

  /* use a reference when you want to update the param directly */
  int* const currCycle = problem->GetView("cycle")->GetNode().as_int32_ptr();

  int numTotalCycles = *(problem->GetView("numTotalCycles")->GetNode().as_int32_ptr());
  int dumpInterval = *(problem->GetView("numCyclesPerDump")->GetNode().as_int32_ptr());

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

