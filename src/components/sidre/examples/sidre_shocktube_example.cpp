/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

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

// Standard library headers
#include <stdio.h>
#include <math.h>
#include <cstdlib>

// SiDRe component headers
#include "sidre/sidre.hpp"

using asctoolkit::sidre::DataBuffer;
using asctoolkit::sidre::DataGroup;
using asctoolkit::sidre::DataStore;
using asctoolkit::sidre::DataView;

using namespace conduit;

#define UPWIND   0
#define DOWNWIND 1

const double gammaa = M_SQRT2;
const double gammaaInverse = M_SQRT1_2;


void CreateScalarIntBufferViewAndSetVal( DataGroup * const grp, const std::string& name, int32 const value )
{
  DataBuffer * const buffer = grp->getDataStore()->createBuffer()
                              ->declare(DataType::int32())
                              ->allocate();


  DataView * const view = grp->createView(name, buffer)->apply(DataType::int32());
  view->setValue(value);

}


void CreateScalarFloatBufferViewAndSetVal( DataGroup * const grp, const std::string& name, float64 const value )
{
  DataBuffer * const buffer = grp->getDataStore()->createBuffer()
                              ->declare(DataType::float64())
                              ->allocate();


  DataView * const view = grp->createView(name,buffer)->apply(DataType::float64());
  view->setValue(value);

}




/**************************************************************************
 * Subroutine:  GetUserInput
 * Purpose   :  Ask for control and output information
 *************************************************************************/

void GetUserInput(DataGroup * const prob)
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
    CreateScalarIntBufferViewAndSetVal( prob, "numElems", numElems );

    CreateScalarIntBufferViewAndSetVal( prob, "numFaces", numFaces );

  }

  /********************/
  /* Get physics info */
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

    CreateScalarFloatBufferViewAndSetVal( prob, "pressureRatio", pratio );

    CreateScalarFloatBufferViewAndSetVal( prob, "densityRatio", dratio );


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

    CreateScalarIntBufferViewAndSetVal( prob, "numUltraDumps", numUltraDumps );

    CreateScalarIntBufferViewAndSetVal( prob, "numCyclesPerDump", numCyclesPerDump );

    CreateScalarIntBufferViewAndSetVal( prob, "numTotalCycles", numUltraDumps * numCyclesPerDump );


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
  int32 const numElems = prob->getView("numElems")->getValue();
  int32 const numFaces = prob->getView("numFaces")->getValue();

  //int inflow[1];
  //int outflow[1];

  /* create element and face classes */

  DataGroup * const elem = prob->createGroup("elem");
  DataGroup * const face = prob->createGroup("face");

  /* set up some important views */

  //inflow[0] = 0; /* identify inflow elements */
//  elem->viewCreate("inflow", new IndexSet(1, inflow));
  elem->createGroup("inflow");

  //outflow[0] = numElems - 1; /* identify outflow elements */
//  elem->viewCreate("outflow", new IndexSet(1, outflow));
  elem->createGroup("outflow");



  int32 numTubeElems = numElems - 2;
  DataGroup * const tube = elem->createGroup("tube");//->SetDataShape(DataStoreNS::DataShape(numTubeElems));

#if 0
  DataBuffer * const mapToElemsBuffer = elem->getDataStore()->createBuffer()
                                        ->declare(DataType::int32(numTubeElems))
                                        ->allocate();


  DataView * const mapToElemsView = tube->createView("mapToElems", mapToElemsBuffer)
                                    ->apply(DataType::int32(numTubeElems));
  int32 * const mapToElems = mapToElemsView->getValue();
#else
  int32 * const mapToElems = tube->createViewAndBuffer("mapToElems", DataType::int32(numTubeElems))
                             ->allocate()
                             ->getValue();
#endif


  for ( int k = 0u ; k < numTubeElems ; ++k)
  {
    mapToElems[k] = k + 1;
  }

  /* Set up some important data relations */

  /* Each face connects to two elements */

#if 0
  DataBuffer * const faceToElemBuffer = face->getDataStore()->createBuffer()
                                        ->declare(DataType::int32(2*numFaces,0,8))
                                        ->allocate();

  int32 * const faceToElem = face->createView("faceToElem", faceToElemBuffer)
                             ->apply(DataType::int32(2*numFaces,0,8))
                             ->getValue();
#else
  int32 * const faceToElem = face->createViewAndBuffer("faceToElem", DataType::int32(2*numFaces))
                             ->allocate()
                             ->getValue();
#endif

  for (i = 0 ; i < numFaces ; ++i)
  {
    faceToElem[i * 2 + UPWIND] = i;
    faceToElem[i * 2 + DOWNWIND] = i + 1;
  }

  /* Each element connects to two faces */ //
//  Relation &elemToFace = *tube->relationCreate("elemToFace", 2);
#if 0
  DataBuffer * const elemToFaceBuffer = tube->getDataStore()->createBuffer()
                                        ->declare(DataType::int32(2*numElems,0,8))
                                        ->allocate();


  int32 * elemToFace = tube->createView("elemToFace", elemToFaceBuffer)
                       ->apply(DataType::int32(2*numElems,0,8))
                       ->getValue();
#else
  int32 * elemToFace = tube->createViewAndBuffer("elemToFace", DataType::int32(2*numElems))
                       ->allocate()
                       ->getValue();
#endif

  for (i = 0 ; i < numElems ; ++i)
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
  DataGroup * const elem = (prob->getGroup("elem"));
  DataGroup * const face = (prob->getGroup("face"));

  /* Create element centered quantities */


  DataBuffer * buffer = ATK_NULLPTR;

  int32 const numElems = prob->getView("numElems")->getValue();
  int32 const numFaces = prob->getView("numFaces")->getValue();

#if 0
  buffer = elem->getDataStore()->createBuffer()
           ->declare(DataType::float64(numElems))
           ->allocate();


  float64 * const mass = elem->createView("mass", buffer)
                         ->apply(DataType::float64(numElems))
                         ->getValue();
#else
  float64 * const mass = elem->createViewAndBuffer("mass", DataType::float64(numElems))
                         ->allocate()
                         ->getValue();
#endif


#if 0
  buffer = elem->getDataStore()->createBuffer()
           ->declare(DataType::float64(numElems))
           ->allocate();

  float64 * const momentum = elem->createView("momentum", buffer)
                             ->apply(DataType::float64(numElems))
                             ->getValue();
#else
  float64 * const momentum = elem->createViewAndBuffer("momentum", DataType::float64(numElems))
                             ->allocate()
                             ->getValue();
#endif

#if 0
  buffer = elem->getDataStore()->createBuffer()
           ->declare(DataType::float64(numElems))
           ->allocate();

  float64 * const energy = elem->createView("energy", buffer)
                           ->apply(DataType::float64(numElems))
                           ->getValue();
#else
  float64 * const energy = elem->createViewAndBuffer("energy", DataType::float64(numElems))
                           ->allocate()
                           ->getValue();
#endif

#if 0
  buffer = elem->getDataStore()->createBuffer()
           ->declare(DataType::float64(numElems))
           ->allocate();

  float64 * const pressure = elem->createView("pressure", buffer)
                             ->apply(DataType::float64(numElems))
                             ->getValue();
#else
  float64 * const pressure = elem->createViewAndBuffer("pressure", DataType::float64(numElems))
                             ->allocate()
                             ->getValue();
#endif

  /* Create face centered quantities */

  buffer = face->getDataStore()->createBuffer()
           ->declare(DataType::float64(numFaces))
           ->allocate();

  face->createView("F0", buffer)
  ->apply(DataType::float64(numFaces));

  buffer = face->getDataStore()->createBuffer()
           ->declare(DataType::float64(numFaces))
           ->allocate();

  face->createView("F1", buffer)
  ->apply(DataType::float64(numFaces));

  buffer = face->getDataStore()->createBuffer()
           ->declare(DataType::float64(numFaces))
           ->allocate();

  face->createView("F2", buffer)
  ->apply(DataType::float64(numFaces));


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
  for (i = startTube ; i < midTube ; ++i)
  {
    mass[i] = massInitial;
    momentum[i] = momentumInitial;
    pressure[i] = pressureInitial;
    energy[i] = energyInitial;
  }

  /* adjust parameters for low pressure portion of tube */
  double dratio = prob->getView("densityRatio")->getValue();
  double pratio = prob->getView("pressureRatio")->getValue();

  massInitial *= dratio;
  pressureInitial *= pratio;
  energyInitial = pressureInitial / (gammaa - 1.0);

  for (i = midTube ; i < endTube ; ++i)
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
  double dx = prob->getView("dx")->getValue();
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

void ComputeFaceInfo(DataGroup * const prob)
{
  int i;
  DataGroup * const face = prob->getGroup("face");
//  Relation &faceToElem = *face->relation("faceToElem");
  int32 const * const faceToElem = face->getView("faceToElem")->getValue();

  float64 * const F0 = face->getView("F0")->getValue();
  float64 * const F1 = face->getView("F1")->getValue();
  float64 * const F2 = face->getView("F2")->getValue();
  int numFaces = face->getView("F0")->getNumberOfElements();

  DataGroup * const elem = prob->getGroup("elem");
  float64 * const mass = elem->getView("mass")->getValue();
  float64 * const momentum = elem->getView("momentum")->getValue();
  float64 * const energy = elem->getView("energy")->getValue();

  for (i = 0 ; i < numFaces ; ++i)
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

void UpdateElemInfo(DataGroup * const prob)
{
  int i;

  /* get the element quantities we want to update */
  DataGroup * const elem = prob->getGroup("elem");
  float64 * const mass = elem->getView("mass")->getValue();
  float64 * const momentum = elem->getView("momentum")->getValue();
  float64 * const energy = elem->getView("energy")->getValue();
  float64 * const pressure = elem->getView("pressure")->getValue();

  /* focus on just the elements within the shock tube */
  DataGroup * const tube = elem->getGroup("tube");
  int32 * const elemToFace = tube->getView("elemToFace")->getValue();

//  Relation &elemToFace = *tube->relation("elemToFace");
  int numTubeElems = tube->getView("mapToElems")->getNumberOfElements();

//  int *is = tube->map();
  int32 * const is = tube->getView("mapToElems")->getValue();

  /* The element update is calculated as the flux between faces */
  DataGroup * const face = prob->getGroup("face");
  float64 * const F0 = face->getView("F0")->getValue();
  float64 * const F1 = face->getView("F1")->getValue();
  float64 * const F2 = face->getView("F2")->getValue();

  double const dx = prob->getView("dx")->getValue();
  double const dt = prob->getView("dt")->getValue();
  float64& time = *prob->getView("time")->getValue<float64 *>();

  for (i = 0 ; i < numTubeElems ; ++i)
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
  FILE * fp;
  char fname[100];
  char * tail;

//   VHashTraverse_t content ;

  strcpy(fname, "problem" );

  /* Skip past the junk */
  for (tail=fname ; isalpha(*tail) ; ++tail)
    ;

  sprintf(tail, "_%04d.ult", prob->getView("cycle")->getValue<int>());

  if ((fp = fopen(fname, "w")) == NULL)
  {
    printf("Could not open file %s. Aborting.\n", fname);
    exit (-1);
  }

  fprintf(fp, "# Ultra Plot File\n");
  fprintf(fp, "# Problem: %s\n", "problem" );


  for(size_t i=0 ; i<prob->getNumViews() ; i++)
  {
    DataView * const view = prob->getView(i);
    const int length = view->getNumberOfElements();
    const std::string& name = view->getName();
    if( length <= 1 )
    {
      if( view->getTypeID() == CONDUIT_INT32_T )
      {
        fprintf(fp, "# %s = %d\n",
                name.c_str(),
                view->getValue<int>());
      }
      else if( view->getTypeID() == CONDUIT_FLOAT64_T )
      {
        fprintf(fp, "# %s = %f\n",
                name.c_str(),
                view->getValue<double>());
      }
    }
  }


  DataGroup * const elem = prob->getGroup("elem");

  for(size_t i=0 ; i<elem->getNumViews() ; i++)
  {
    DataView * const view = elem->getView(i);
    const int length = view->getNumberOfElements();
    const std::string& name = view->getName();
    fprintf(fp, "# %s\n", name.c_str() );

    if( view->getTypeID() == CONDUIT_INT32_T )
    {
      int32 const * const data = view->getValue();
      for ( int j=0 ; j<length ; ++j)
      {
        fprintf(fp, "%f %f\n", (double) j, (double) data[j]);
      }
      fprintf(fp, "\n");
    }
    else if( view->getTypeID() == CONDUIT_FLOAT64_T )
    {
      float64 const * const data = view->getValue();
      for ( int j=0 ; j<length ; ++j)
      {
        fprintf(fp, "%f %f\n", (double) j, (double) data[j]);
      }
      fprintf(fp, "\n");
    }
  }


  fclose(fp);
#endif
  return;
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

  DataStore DATASTORE;
  DataStore * const dataStore = &DATASTORE;
  DataGroup * const rootGroup = dataStore->getRoot();
  DataGroup * const prob = rootGroup->createGroup("problem");

  GetUserInput(prob);
  CreateShockTubeMesh(prob);
  InitializeShockTube(prob);

  /* use a reference when you want to update the param directly */
  int * const currCycle = prob->getView("cycle")->getValue();

  int numTotalCycles = prob->getView("numTotalCycles")->getValue();
  int dumpInterval = prob->getView("numCyclesPerDump")->getValue();

  for (*currCycle = 0 ; *currCycle < numTotalCycles ; ++(*currCycle) )
  {
    /* dump the ultra file, based on the user chosen attribute mask */
    if ( (*currCycle) % dumpInterval == 0)
    {
      DumpUltra(prob);
    }

    ComputeFaceInfo(prob);
    UpdateElemInfo(prob);
  }

  DumpUltra(prob); /* One last dump */

  return 0;
}
