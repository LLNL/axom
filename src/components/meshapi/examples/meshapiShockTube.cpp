/**
 * \file meshapiShockTube.cc
 *
 * \brief 1D shock tube, split flux Euler equations
 *
 * \author J. Keasler (original)
 * \author K. Weiss (modified to use the ASC Toolkit Mesh API)
 *
 * \details
 * \verbatim
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
 *     All quantities are non-dimensionalized.
 *
 *     @Q   @F    @Q   @F @Q
 *     -- + -- =  -- + -- -- = 0
 *     @t   @x    @t   @Q @x
 * \endverbatim
 */


#include <cmath>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>


#include "common/Utilities.hpp"
#include "meshapi/FieldRegistry.hpp"
#include "meshapi/OrderedSet.hpp"
#include "meshapi/StaticConstantRelation.hpp"


namespace asctoolkit {
namespace meshapi {
namespace examples {
namespace shocktube {

    asctoolkit::meshapi::MeshIndexType const UPWIND   = 0;
    asctoolkit::meshapi::MeshIndexType const DOWNWIND = 1;

    const double gammaa = M_SQRT2 ;
    const double gammaaInverse = M_SQRT1_2 ;

    const int INIT_NUM_ELEMS = 100;
    const int INIT_NUM_ULTRA_DUMPS = 5;
    const int INIT_NUM_ULTRA_CYCLES_PER_DUMP = 200;

    const double INIT_P_RATIO = 0.5;
    const double INIT_D_RATIO = 0.5;

    /**
     * \brief Simple representation of the mesh for this 1D example
     *
     * \details Mesh contains a set of elements and a set of faces between elements
     *         as well as the relations from elements to faces and vice versa.
     *
     * \note The mesh is currently missing a few subsets on the element set
     *       Specifically, we want the 'tube' to skip the first and last element
     *       ('inflow' and 'outflow' in the code below)
     *
     * \note We are also missing an implicit constant grid relation.
     * \note For current implementation with explicit static (constant) relations.
     *       We are missing a nice way to set the relation elements.
     *       It should not have to be done explicitly in user code
     *       Idea: We could have a relationInverter function that takes a relation from sets A to B
     *             and generates a relation from set B to set A with all the arrows reversed.
     */
    class ShockTubeMesh
    {
    public:
        // types for sets
        typedef asctoolkit::meshapi::OrderedSet ElemSet;
        typedef asctoolkit::meshapi::OrderedSet FaceSet;

        // types for relations
        typedef asctoolkit::meshapi::StaticConstantRelation ElemToFaceRelation;
        typedef asctoolkit::meshapi::StaticConstantRelation FaceToElemRelation;

        // other types
        typedef asctoolkit::meshapi::Set::SetIndex IndexType;
        typedef asctoolkit::meshapi::Set::SizeType SizeType;

    public:
        ElemSet elems;
        FaceSet faces;

        FaceToElemRelation relationFaceElem;
        ElemToFaceRelation relationElemFace;

    };


    // Define the explicit instances for int and double
    FieldRegistry<int>    intsRegistry;
    FieldRegistry<double> realsRegistry;

    typedef FieldRegistry<int>::MapType    IntField;
    typedef FieldRegistry<double>::MapType RealField;


/**************************************************************************
 * Subroutine:  GetUserInput
 * Purpose   :  Ask for control and output information
 *************************************************************************/

void GetUserInput()
{
   /**********************************/
   /* Get mesh info, and create mesh */
   /**********************************/
   {
      int numElems;
      int numFaces;

      std::cout<< "How many zones for the 1D shock tube? ";
      numElems = INIT_NUM_ELEMS;
      std::cout << numElems << std::endl;

      // add an inflow and outflow zone
      numElems += 2;
      numFaces = numElems - 1;

      intsRegistry.addScalar("numElems", numElems);
      intsRegistry.addScalar("numFaces", numFaces) ;
   }

   /********************/
   /* Get pyhsics info */
   /********************/
   {
      double pratio = -1.0;
      double dratio = -1.0 ;

      while (pratio < 0.0 || pratio > 1.0) {
         std::cout<< "What pressure ratio would you like (0 <= x <= 1)? ";
         pratio = INIT_P_RATIO;
         std::cout<< pratio <<std::endl;
      }

      while (dratio < 0.0 || dratio > 1.0) {
         std::cout << "What density ratio would you like (0 <= x <= 1)? ";
         dratio = INIT_D_RATIO;
         std::cout << dratio <<std::endl;
      }

      realsRegistry.addScalar("pressureRatio", pratio) ;
      realsRegistry.addScalar("densityRatio", dratio) ;
   }

   /********************/
   /* Get output  info */
   /********************/
   {
      int numUltraDumps;
      int numCyclesPerDump ;

      std::cout << "How many Ultra dumps would you like? ";
      numUltraDumps = INIT_NUM_ULTRA_DUMPS;
      std::cout << numUltraDumps <<std::endl;

      std::cout << "How many cycles per Ultra dump would you like? ";
      numCyclesPerDump = INIT_NUM_ULTRA_CYCLES_PER_DUMP;
      std::cout << numCyclesPerDump <<std::endl;

      int numTotalCycles = numUltraDumps*numCyclesPerDump;
      std::cout << "\nSimulation will run for " << numTotalCycles <<" cycles."<<std::endl;

      intsRegistry.addScalar("numUltraDumps", numUltraDumps) ;
      intsRegistry.addScalar("numCyclesPerDump", numCyclesPerDump) ;
      intsRegistry.addScalar("numTotalCycles", numTotalCycles) ;
   }

   return ;
}

/**
 * \brief Build an empty mesh for the shock tube
 * \details
 * \verbatim
 *      Gaps between elements are faces
 *                     |
 *      -------------------------------
 *      |   |   |             |   |   |
 *
 *   ### ### ### ###       ### ### ### ###
 *   ### ### ### ###  ...  ### ### ### ###  <--- 1D Shock tube model
 *   ### ### ### ###       ### ### ### ###
 *
 *    |  |                           |  |
 *    |  -----------------------------  |
 *   Inflow           |               Outflow
 *   Element      Tube Elements       Element
 *
 * \endverbatim
 */
void CreateShockTubeMesh(ShockTubeMesh *mesh)
{
   // create element and face sets
   mesh->elems = ShockTubeMesh::ElemSet( intsRegistry.getScalar("numElems") );
   mesh->faces = ShockTubeMesh::FaceSet( intsRegistry.getScalar("numFaces") );

   // TODO: Need to define subsets for inflow, outflow and tube

   /*
       int inflow[1]  ;
       int outflow[1] ;
       inflow[0] = 0 ;              // identify inflow elements
       elem->viewCreate("inflow", new IndexSet(1, inflow)) ;

       outflow[0] = numElems - 1 ;  // identify outflow elements
       elem->viewCreate("outflow", new IndexSet(1, outflow)) ;

       // identify shock tube elements - set up basic map
       View *tube = elem->viewCreate("tube", new IndexSet(numElems-2)) ;

       // Shift IndexSet indices to identify shocktube elements
       // (shocktube element numbers are one through numElems-1 inclusive)
       tube->indexSet()->shift(1) ;
   */

   // ------------ Set up relations

   // TODO: Need to define DynamicConstantRelation -- which will allow modifying the relation elements
   // TODO: Need to define ImplicitConstantRelation -- since the relations are actually implicit
   //       -- no storage should be necessary for regular grid neighbors
   // For now, we will have to do this explicitly...
   const int STRIDE = 2;
   typedef std::vector<ShockTubeMesh::IndexType> IndexVec;

   /// Setup the face -> elem relation
   IndexVec relVec( STRIDE * mesh->faces.size());
   IndexVec::iterator relIt = relVec.begin();
   for(ShockTubeMesh::IndexType idx =0; idx < static_cast<ShockTubeMesh::IndexType>(mesh->faces.size()); ++idx)
   {
       *relIt++ = idx;
       *relIt++ = idx+1;
   }

   mesh->relationFaceElem = ShockTubeMesh::FaceToElemRelation(&mesh->faces, &mesh->elems);
   mesh->relationFaceElem.setRelation(relVec, STRIDE);
   ATK_ASSERT(mesh->relationFaceElem.isValid());


   /// Setup the elem -> face relation
   // Note: This relation needs to be on the TUBE subset of the elements (i.e. skip the first and last elt
   //       It is currently on all elements, which necessitates a workaround below.
   //       And the first and last elements of the relation are set to 0 which is wrong.
   //       Question -- how are we planning to handle indexes that are out or range (accidentally)?
   //                   how are we planning to handle indexes that are intentionally out of range
   //                   (e.g. to indicate a problem, or a missing element etc..)?
   ShockTubeMesh::SizeType elemSize = mesh->elems.size();
   relVec.clear();
   relVec.resize( STRIDE * elemSize);
   relIt = relVec.begin();
   for(ShockTubeMesh::IndexType idx =0; idx < static_cast<ShockTubeMesh::IndexType>(elemSize); ++idx)
   {
       // HACK -- these indexes should be over the tube subset -- skipping the first and last elem
       *relIt++ = (idx == 0)? 0 : idx-1;
       *relIt++ = (idx == static_cast<ShockTubeMesh::IndexType>(elemSize-1 ))? 0 : idx;
   }

   mesh->relationElemFace = ShockTubeMesh::ElemToFaceRelation(&mesh->elems, &mesh->faces);
   mesh->relationElemFace.setRelation(relVec, STRIDE);
   ATK_ASSERT(mesh->relationElemFace.isValid());

}


/**************************************************************************
 * Subroutine:  InitializeShockTube
 * Purpose   :  Populate the mesh with values
 *************************************************************************/
void InitializeShockTube(ShockTubeMesh const& mesh)
{
    typedef ShockTubeMesh::IndexType IndexType;

   // TODO: Define and use mesh API maps over sets for these
   // Note -- the extra code and allocations here will be unnecessary once maps on sets are defined

   // Create element centered fields
   RealField& mass      = realsRegistry.addField("mass", &mesh.elems);
   RealField& momentum  = realsRegistry.addField("momentum", &mesh.elems);
   RealField& energy    = realsRegistry.addField("energy", &mesh.elems);
   RealField& pressure  = realsRegistry.addField("pressure", &mesh.elems );

   // Create face centered fields
   // mv, mv^2+P, and v(E+P)
   realsRegistry.addField("F0", &mesh.faces );
   realsRegistry.addField("F1", &mesh.faces );
   realsRegistry.addField("F2", &mesh.faces );

   // Fill left half with high pressure, right half with low pressure
   IndexType startTube  = 0 ;
   IndexType endTube    = mesh.elems.size();
   IndexType midTube    = endTube / 2 ;

   // Non-dimensionalized reference values
   double massInitial       = 1.0 ;
   double momentumInitial   = 0.0 ;
   double pressureInitial   = gammaaInverse ;
   double energyInitial     = pressureInitial/(gammaa-1.0) ;

   // Initialize zonal quantities
   for (IndexType idx=startTube; idx<midTube; ++idx)
   {
      mass[idx]     = massInitial ;
      momentum[idx] = momentumInitial ;
      pressure[idx] = pressureInitial ;
      energy[idx]   = energyInitial ;
   }

   // adjust parameters for low pressure portion of tube
   double dratio = realsRegistry.getScalar("densityRatio") ;
   double pratio = realsRegistry.getScalar("pressureRatio") ;

   massInitial      *= dratio ;
   pressureInitial  *= pratio ;
   energyInitial     = pressureInitial/(gammaa - 1.0) ;

   for (IndexType idx=midTube; idx<endTube; ++idx)
   {
      mass[idx]     = massInitial ;
      momentum[idx] = momentumInitial ;
      pressure[idx] = pressureInitial ;
      energy[idx]   = energyInitial ;
   }

   // Create needed time info
   realsRegistry.addScalar("time", 0.0) ;
   intsRegistry.addScalar("cycle", 0) ;

   double dx = 1.0 / static_cast<double>(endTube);
   realsRegistry.addScalar("dx", dx);
   realsRegistry.addScalar("dt", 0.4*dx) ;

}

/**
 * \function ComputeFaceInfo
 * \brief Compute F quantities at faces.
 *
 * \details
 *
 * \verbatim
 *  @F   @F0   @F1   @F2
 *  -- = --- + --- + ---  
 *  @x   @x    @x    @x
 *  \endverbatime
 *
 *  Calculate F0, F1 and F2 at the face centers.
 *
 */
void ComputeFaceInfo(ShockTubeMesh const& mesh)
{
   typedef ShockTubeMesh::IndexType IndexType;

   // Face fields
   RealField & F0 = realsRegistry.getField("F0") ;
   RealField & F1 = realsRegistry.getField("F1") ;
   RealField & F2 = realsRegistry.getField("F2") ;

   // Element fields
   RealField const& mass     = realsRegistry.getField("mass") ;
   RealField const& momentum = realsRegistry.getField("momentum") ;
   RealField const& energy   = realsRegistry.getField("energy") ;

   // Update face data using element data using the face->elem relation
   for (ShockTubeMesh::IndexType fIdx=0; fIdx< static_cast<ShockTubeMesh::IndexType>(mesh.faces.size() ); ++fIdx)
   {
      // each face has an upwind and downwind element.
      IndexType upWind   = mesh.relationFaceElem[fIdx][UPWIND] ;   // upwind element
      IndexType downWind = mesh.relationFaceElem[fIdx][DOWNWIND] ; // downwind element

      // calculate face centered quantities as avg of element centered ones
      double massf      = 0.5 * (mass[upWind]     + mass[downWind] ) ;
      double momentumf  = 0.5 * (momentum[upWind] + momentum[downWind] ) ;
      double energyf    = 0.5 * (energy[upWind]   + energy[downWind] ) ;
      double pressuref  = (gammaa - 1.0) * (energyf - 0.5*momentumf*momentumf/massf) ;
      double c = sqrt(gammaa*pressuref/massf) ;
      double v = momentumf/massf ;

      double ev ; 
      double cLocal ;


      // Now that we have the wave speeds, we might want to
      // look for the max wave speed here, and update dt
      // appropriately right before leaving this function.

      // OK, calculate face quantities

      F0[fIdx] = F1[fIdx] = F2[fIdx] = 0.0 ;

      IndexType contributor = ((v >= 0.0) ? upWind : downWind) ;
      massf     = mass[contributor] ;
      momentumf = momentum[contributor] ;
      energyf   = energy[contributor] ;
      pressuref = energyf - 0.5*momentumf*momentumf/massf ;
      ev        = v*(gammaa - 1.0) ;

      F0[fIdx] += ev*massf ;
      F1[fIdx] += ev*momentumf ;
      F2[fIdx] += ev*(energyf - pressuref) ;

      contributor = ((v + c >= 0.0) ? upWind : downWind) ;
      massf     = mass[contributor] ;
      momentumf = momentum[contributor] ;
      energyf   = energy[contributor] ;
      pressuref = (gammaa - 1.0)*(energyf - 0.5*momentumf*momentumf/massf) ;
      ev        = 0.5*(v + c) ;
      cLocal    = sqrt(gammaa*pressuref/massf) ;

      F0[fIdx] += ev*massf ;
      F1[fIdx] += ev*(momentumf + massf*cLocal) ;
      F2[fIdx] += ev*(energyf + pressuref + momentumf*cLocal) ;

      contributor = ((v - c >= 0.0) ? upWind : downWind) ;
      massf     = mass[contributor] ;
      momentumf = momentum[contributor] ;
      energyf   = energy[contributor] ;
      pressuref = (gammaa - 1.0)*(energyf - 0.5*momentumf*momentumf/massf) ;
      ev        = 0.5*(v - c) ;
      cLocal    = sqrt(gammaa*pressuref/massf) ;

      F0[fIdx] += ev*massf ;
      F1[fIdx] += ev*(momentumf - massf*cLocal) ;
      F2[fIdx] += ev*(energyf + pressuref - momentumf*cLocal) ;
   }
}

/**************************************************************************
 * Subroutine:  UpdateElemInfo
 * Purpose   :  Q(elem) = Q(elem) + deltaQ(elem)
 *
 *  deltaQ(elem) = - (F(downWindFace) - F(upWindFace)) * dt / dx ;
 *
 *************************************************************************/

void UpdateElemInfo(ShockTubeMesh const& mesh)
{
    // get the element quantities we want to update
    RealField & mass     = realsRegistry.getField("mass") ;
    RealField & momentum = realsRegistry.getField("momentum") ;
    RealField & energy   = realsRegistry.getField("energy") ;
    RealField & pressure = realsRegistry.getField("pressure") ;

   // The element update is calculated as the flux between faces
   RealField const& F0 = realsRegistry.getField("F0") ;
   RealField const& F1 = realsRegistry.getField("F1") ;
   RealField const& F2 = realsRegistry.getField("F2") ;

   double   dx   = realsRegistry.getScalar("dx") ;
   double   dt   = realsRegistry.getScalar("dt") ;
   double & time = realsRegistry.getScalar("time") ;


   /// Update the element fields based on the face data using the elem->face relation
   ShockTubeMesh::ElemSet::iterator_pair elemItPair = mesh.elems.range();
   // HACK: We must update the iterator ranges to match the tube.
   // TODO: Switch to relation on tube element set once Subsets are defined.
   elemItPair.first++;    elemItPair.second--;

   for (; elemItPair.first < elemItPair.second; ++elemItPair.first)
   {
      ShockTubeMesh::IndexType elemIdx = *elemItPair.first;

      // Each element inside the tube has an upwind and downwind face
      ShockTubeMesh::IndexType upWind   = mesh.relationElemFace[elemIdx][UPWIND] ;      // upwind face
      ShockTubeMesh::IndexType downWind = mesh.relationElemFace[elemIdx][DOWNWIND] ;    // downwind face

      mass[elemIdx]     -= gammaaInverse*(F0[downWind] - F0[upWind])*dt/dx ;
      momentum[elemIdx] -= gammaaInverse*(F1[downWind] - F1[upWind])*dt/dx ;
      energy[elemIdx]   -= gammaaInverse*(F2[downWind] - F2[upWind])*dt/dx ;
      pressure[elemIdx]  = (gammaa - 1.0) * (energy[elemIdx] - 0.5 * momentum[elemIdx]*momentum[elemIdx]/mass[elemIdx]) ;
   }

   // update the time
   time += dt ;
}

void dumpData(ShockTubeMesh const& mesh)
{
    RealField const& mass     = realsRegistry.getField("mass") ;
    RealField const& momentum = realsRegistry.getField("momentum") ;
    RealField const& energy   = realsRegistry.getField("energy") ;
    RealField const& pressure = realsRegistry.getField("pressure") ;

    static const ShockTubeMesh::ElemSet::SizeType MAX_ELEM_DUMP = 10;
    const int maxDump = std::min(mesh.elems.size(), MAX_ELEM_DUMP);
    const int rmaxDump = std::min(MAX_ELEM_DUMP, mesh.elems.size() - maxDump);

    // TODO: The following is currently grabbing the raw data from the Map and spitting out at most MAX_ELEM_DUMP elements
    // I would like to create a subset with a stride to only print every n_th element
    std::cout<<"\n\t\tElem idx: ";
    std::copy(mesh.elems.begin(), mesh.elems.begin()+maxDump, std::ostream_iterator<ShockTubeMesh::IndexType>(std::cout, "\t"));
    std::cout<<"...\t";
    std::copy(mesh.elems.end()-rmaxDump, mesh.elems.end(), std::ostream_iterator<ShockTubeMesh::IndexType>(std::cout, "\t"));

    std::cout<<"\n\t\tMass    : " << std::setprecision(3) ;
    std::copy(mass.data().begin(), mass.data().begin()+maxDump, std::ostream_iterator<double>(std::cout, "\t"));
    std::cout<<"...\t";
    std::copy(mass.data().end()-rmaxDump, mass.data().end(), std::ostream_iterator<double>(std::cout, "\t"));

    std::cout<<"\n\t\tMomentum: ";
    std::copy(momentum.data().begin(), momentum.data().begin()+maxDump, std::ostream_iterator<double>(std::cout, "\t"));
    std::cout<<"...\t";
    std::copy(momentum.data().end()-rmaxDump, momentum.data().end(), std::ostream_iterator<double>(std::cout, "\t"));

    std::cout<<"\n\t\tEnergy  : ";
    std::copy(energy.data().begin(), energy.data().begin()+maxDump, std::ostream_iterator<double>(std::cout, "\t"));
    std::cout<<"...\t";
    std::copy(energy.data().end()-rmaxDump, energy.data().end(), std::ostream_iterator<double>(std::cout, "\t"));

    std::cout<<"\n\t\tPressure: ";
    std::copy(pressure.data().begin(), pressure.data().begin()+maxDump, std::ostream_iterator<double>(std::cout, "\t"));
    std::cout<<"...\t";
    std::copy(pressure.data().end()-rmaxDump, pressure.data().end(), std::ostream_iterator<double>(std::cout, "\t"));

}

} // end namespace shocktube
} // end namespace examples
} // end namespace meshapi
} // end namespace asctoolkit


/**************************************************************************
 * Subroutine:  main
 * Purpose   :  Simulate a 1D Shock Tube using split flux Euler formulation
 *************************************************************************/

int main(void)
{
    using namespace asctoolkit::meshapi::examples::shocktube;

   //extern void DumpUltra(View *prob) ;

   // We should be able to parallelize pretty easily by
   // adding an MPI_Init() here, modifying the setup slightly,
   // adding a communication subroutine in the main loop,
   // and calling MPI_Finalize() at the end of main()

   //View *problem = new View("ShockTube, 1D Split Flux Euler") ;

    ShockTubeMesh mesh;

    GetUserInput();

    CreateShockTubeMesh(&mesh);     // setup sets and relations
    InitializeShockTube(mesh);      // setup fields

    int numTotalCycles = intsRegistry.getScalar("numTotalCycles") ;
    int dumpInterval   = intsRegistry.getScalar("numCyclesPerDump") ;

    // use the & operation when you want to update the param directly
    int& currCycle = intsRegistry.getScalar("cycle") ;
    for (currCycle=0; currCycle<numTotalCycles; ++currCycle)
    {
        if( currCycle % dumpInterval == 0)
        {
            std::cout<< "\n\tStarting cycle " << currCycle << " at time " << realsRegistry.getScalar("time");
            dumpData(mesh);
        }

        ComputeFaceInfo(mesh) ;
        UpdateElemInfo(mesh) ;
    }

    std::cout<< "\n\tFinished cycle " << currCycle << " at time " << realsRegistry.getScalar("time");
    dumpData(mesh);
    std::cout<<"\ndone." << std::endl;

    return 0 ;
}
