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


#include "Vista.h"
#include "View.h"

#include <stdio.h>
#include <math.h>
#define UPWIND   0
#define DOWNWIND 1

const double gammaa = M_SQRT2 ;
const double gammaaInverse = M_SQRT1_2 ;

/**************************************************************************
 * Subroutine:  GetUserInput
 * Purpose   :  Ask for control and output information
 *************************************************************************/

void GetUserInput(View *problem)
{
   /**********************************/
   /* Get mesh info, and create mesh */
   /**********************************/
   {
      int numElems, numFaces ;

      printf("How many zones for the 1D shock tube? ") ;
      scanf("%d", &numElems) ;

      numElems += 2 ;            /* add an inflow and outflow zone */
      numFaces = numElems - 1 ;

      problem->paramAddInt("numElems", numElems) ;
      problem->paramAddInt("numFaces", numFaces) ;
   }

   /********************/
   /* Get pyhsics info */
   /********************/
   {
      double pratio = -1.0 , dratio = -1.0 ;

#if 0
      printf("What cfl number would you like to use? ") ;
      scanf("%lf", &cfl) ;
#endif

      while (pratio < 0.0 || pratio > 1.0) {
         printf("What pressure ratio would you like (0 <= x <= 1)? ") ;
         scanf("%lf", &pratio) ;
      }

      while (dratio < 0.0 || dratio > 1.0) {
         printf("What density ratio would you like (0 <= x <= 1)? ") ;
         scanf("%lf", &dratio) ;
      }

      problem->paramAddReal("pressureRatio", pratio) ;
      problem->paramAddReal("densityRatio", dratio) ;
   }

   /********************/
   /* Get output  info */
   /********************/
   {
      int numUltraDumps, numCyclesPerDump ;

      printf("How many Ultra dumps would you like? ") ;
      scanf("%d", &numUltraDumps) ;

      printf("How many cycles per Ultra dump would you like? ") ;
      scanf("%d", &numCyclesPerDump) ;

      problem->paramAddInt("numUltraDumps", numUltraDumps) ;
      problem->paramAddInt("numCyclesPerDump", numCyclesPerDump) ;
      problem->paramAddInt("numTotalCycles",
                            numUltraDumps*numCyclesPerDump) ;
   }

   return ;
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

void CreateShockTubeMesh(View *prob)
{
   int i ;
   int numElems = prob->paramInt("numElems") ;
   int numFaces = prob->paramInt("numFaces") ;
   int inflow[1]  ;
   int outflow[1] ;

   /* create element and face classes */

   View *elem = prob->viewCreate("elem", new IndexSet(numElems)) ;
   View *face = prob->viewCreate("face", new IndexSet(numFaces)) ;

   /* set up some important views */

   inflow[0] = 0 ;              /* identify inflow elements */
   elem->viewCreate("inflow", new IndexSet(1, inflow)) ;

   outflow[0] = numElems - 1 ;  /* identify outflow elements */
   elem->viewCreate("outflow", new IndexSet(1, outflow)) ;

   /* identify shock tube elements - set up basic map */
   View *tube = elem->viewCreate("tube", new IndexSet(numElems-2)) ;

   /* Shift IndexSet indicies to identify shocktube elements */
   /* (shocktube element numbers are one through numElems-1 inclusive) */
   tube->indexSet()->shift(1) ;

   /* Set up some important data relations */

   /* Each face connects to two elements */
   Relation &faceToElem = *face->relationCreate("faceToElem", 2) ;
   for (i=0; i<faceToElem.length(); ++i)
   {
      faceToElem[i][UPWIND] = i ;
      faceToElem[i][DOWNWIND] = i+1 ;
   }

   /* Each element connects to two faces */
   Relation &elemToFace = *tube->relationCreate("elemToFace", 2) ;
   for (i=0; i<elemToFace.length(); ++i)
   {
      elemToFace[i][UPWIND] = i ;   /* same map as above by coincidence */
      elemToFace[i][DOWNWIND] = i+1 ;
   }

   return ;
}

/**************************************************************************
 * Subroutine:  InitializeShockTube
 * Purpose   :  Populate the mesh with values
 *************************************************************************/

void InitializeShockTube(View *prob)
{
   int i ;

   /* These were created in GetUserInput() */
   View *elem = prob->view("elem") ;
   View *face = prob->view("face") ;

   /* Create element centered quantities */
   double *mass =     elem->fieldCreateReal("mass")->Real() ;
   double *momentum = elem->fieldCreateReal("momentum")->Real() ;
   double *energy =   elem->fieldCreateReal("energy")->Real() ;
   double *pressure = elem->fieldCreateReal("pressure")->Real() ;

   /* Create face centered quantities */
   face->fieldCreateReal("F", 3) ;  /* mv, mv^2+P, and v(E+P) */

   /* Fill left half with high pressure, right half with low pressure */
   int startTube = 0 ;
   int endTube = elem->length() ;
   int midTube = endTube / 2 ;

   /* Non-dimensionalized reference values */
   double massInitial = 1.0 ;
   double momentumInitial = 0.0 ;
   double pressureInitial = gammaaInverse ;
   double energyInitial = pressureInitial/(gammaa-1.0) ;

   /* Initialize zonal quantities*/
   for (i=startTube; i<midTube; ++i)
   {
      mass[i] = massInitial ;
      momentum[i] = momentumInitial ;
      pressure[i] = pressureInitial ;
      energy[i] = energyInitial ;
   }

   /* adjust parameters for low pressure portion of tube */
   double dratio = prob->paramReal("densityRatio") ;
   double pratio = prob->paramReal("pressureRatio") ;

   massInitial *= dratio ;
   pressureInitial *= pratio ;
   energyInitial = pressureInitial/(gammaa - 1.0) ;

   for (i=midTube; i<endTube; ++i)
   {
      mass[i] = massInitial ;
      momentum[i] = momentumInitial ;
      pressure[i] = pressureInitial ;
      energy[i] = energyInitial ;
   }

   /* Create needed time info */
   prob->paramAddReal("time", 0.0) ;
   prob->paramAddInt("cycle", 0) ;
   prob->paramAddReal("dx", (1.0 / ((double) endTube))) ;
   double dx = prob->paramReal("dx") ;
   prob->paramAddReal("dt", 0.4*dx) ;

   return ;
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

void ComputeFaceInfo(View *problem)
{
   int i ;
   View *face = problem->view("face") ;
   Relation &faceToElem = *face->relation("faceToElem") ;
   double *F0 = face->fieldReal("F", 0) ;
   double *F1 = face->fieldReal("F", 1) ;
   double *F2 = face->fieldReal("F", 2) ;
   int numFaces = face->length() ;

   View *elem = problem->view("elem") ;
   double *mass =     elem->fieldReal("mass") ;
   double *momentum = elem->fieldReal("momentum") ;
   double *energy =   elem->fieldReal("energy") ;

   for (i=0; i<numFaces; ++i)
   {
      /* each face has an upwind and downwind element. */
      int upWind   = faceToElem[i][UPWIND] ;   /* upwind element */
      int downWind = faceToElem[i][DOWNWIND] ; /* downwind element */

      /* calculate face centered quantities */
      double massf =     0.5 * (mass[upWind]     + mass[downWind] ) ;
      double momentumf = 0.5 * (momentum[upWind] + momentum[downWind] ) ;
      double energyf =   0.5 * (energy[upWind]   + energy[downWind] ) ;
      double pressuref = (gammaa - 1.0) *
                         (energyf - 0.5*momentumf*momentumf/massf) ;
      double c = sqrt(gammaa*pressuref/massf) ;
      double v = momentumf/massf ;
      double ev ; 
      double cLocal ;
      int contributor ;

      /* Now that we have the wave speeds, we might want to */
      /* look for the max wave speed here, and update dt */
      /* appropriately right before leaving this function. */
      /* ... */

      /* OK, calculate face quantities */

      F0[i] = F1[i] = F2[i] = 0.0 ;

      contributor = ((v >= 0.0) ? upWind : downWind) ;
      massf = mass[contributor] ;
      momentumf = momentum[contributor] ;
      energyf = energy[contributor] ;
      pressuref = energyf - 0.5*momentumf*momentumf/massf ;
      ev = v*(gammaa - 1.0) ;

      F0[i] += ev*massf ;
      F1[i] += ev*momentumf ;
      F2[i] += ev*(energyf - pressuref) ;

      contributor = ((v + c >= 0.0) ? upWind : downWind) ;
      massf = mass[contributor] ;
      momentumf = momentum[contributor] ;
      energyf = energy[contributor] ;
      pressuref = (gammaa - 1.0)*(energyf - 0.5*momentumf*momentumf/massf) ;
      ev = 0.5*(v + c) ;
      cLocal = sqrt(gammaa*pressuref/massf) ;

      F0[i] += ev*massf ;
      F1[i] += ev*(momentumf + massf*cLocal) ;
      F2[i] += ev*(energyf + pressuref + momentumf*cLocal) ;

      contributor = ((v - c >= 0.0) ? upWind : downWind) ;
      massf = mass[contributor] ;
      momentumf = momentum[contributor] ;
      energyf = energy[contributor] ;
      pressuref = (gammaa - 1.0)*(energyf - 0.5*momentumf*momentumf/massf) ;
      ev = 0.5*(v - c) ;
      cLocal = sqrt(gammaa*pressuref/massf) ;

      F0[i] += ev*massf ;
      F1[i] += ev*(momentumf - massf*cLocal) ;
      F2[i] += ev*(energyf + pressuref - momentumf*cLocal) ;
   }

   return ;
}


/**************************************************************************
 * Subroutine:  UpdateElemInfo
 * Purpose   :  Q(elem) = Q(elem) + deltaQ(elem)
 *
 *  deltaQ(elem) = - (F(downWindFace) - F(upWindFace)) * dt / dx ;
 *
 *************************************************************************/

void UpdateElemInfo(View *problem)
{
   int i ;

   /* get the element quantities we want to update */
   View *elem =      problem->view("elem") ;
   double *mass =     elem->fieldReal("mass") ;
   double *momentum = elem->fieldReal("momentum") ;
   double *energy =   elem->fieldReal("energy") ;
   double *pressure = elem->fieldReal("pressure") ;

   /* focus on just the elements within the shock tube */
   View *tube = elem->view("tube") ;
   Relation &elemToFace = *tube->relation("elemToFace") ;
   int numTubeElems = tube->length() ;
   int *is = tube->map() ;

   /* The element update is calculated as the flux between faces */
   View *face = problem->view("face") ;
   double *F0 = face->fieldReal("F", 0) ;
   double *F1 = face->fieldReal("F", 1) ;
   double *F2 = face->fieldReal("F", 2) ;

   double dx = problem->paramReal("dx") ;
   double dt = problem->paramReal("dt") ;
   double &time = problem->paramReal("time") ;

   for (i=0; i<numTubeElems; ++i)
   {
      /* recalculate elements in the shocktube, don't touch inflow/outflow */
      int elemIdx = is[i] ;

      /* each element inside the tube has an upwind and downwind face */
      int upWind = elemToFace[i][UPWIND] ;     /* upwind face */
      int downWind = elemToFace[i][DOWNWIND] ; /* downwind face */

      mass[elemIdx] -= gammaaInverse*(F0[downWind] - F0[upWind])*dt/dx ;
      momentum[elemIdx] -= gammaaInverse*(F1[downWind] - F1[upWind])*dt/dx ;
      energy[elemIdx] -= gammaaInverse*(F2[downWind] - F2[upWind])*dt/dx ;
      pressure[elemIdx] = (gammaa - 1.0) *
                          (energy[elemIdx] - 0.5 *
                           momentum[elemIdx]*momentum[elemIdx]/mass[elemIdx]) ;
   }

   /* update the time */
   time += dt ;

   return ;
}


/**************************************************************************
 * Subroutine:  main
 * Purpose   :  Simulate a 1D Shock Tube using split flux Euler formulation
 *************************************************************************/

int main(void)
{
   extern void DumpUltra(View *prob) ;

   /* We should be able to parallelize pretty easily by */
   /* adding an MPI_Init() here, modifying the setup slightly, */
   /* adding a communication subroutine in the main loop, */
   /* and calling MPI_Finalize() at the end of main() */

   View *problem = new View("ShockTube, 1D Split Flux Euler") ;

   GetUserInput(problem) ;
   CreateShockTubeMesh(problem) ;
   InitializeShockTube(problem) ; 

   /* use the & operation when you want to update the param directly */
   int &currCycle = problem->paramInt("cycle") ;

   int numTotalCycles = problem->paramInt("numTotalCycles") ;
   int dumpInterval = problem->paramInt("numCyclesPerDump") ;

   for (currCycle=0; currCycle<numTotalCycles; ++currCycle)
   {
      /* dump the ultra file, based on the user chosen arrtibute mask */
      if (currCycle % dumpInterval == 0)
         DumpUltra(problem) ;

      ComputeFaceInfo(problem) ;
      UpdateElemInfo(problem) ;
   }

   DumpUltra(problem) ; /* One last dump */

   return 0 ;
}
