/**
 * \file meshapiShockTube.cpp
 *
 * \brief 1D shock tube, split flux Euler equations
 *
 * \author J. Keasler (original)
 * \author K. Weiss (modified to use the ASC Toolkit Mesh API)
 *
 * \details  Developing example to use and demo features of Mesh API on shock tube example over structured 1D mesh.
 *           Tests: Sets and subsets.
 *                  Implicit relations over regular grid (currently implemented as explicit static constant relations) 5/2015
 *                  Fields/maps over the data -- and access to sidre/local datastore.
 *
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
#include "meshapi/RangeSet.hpp"
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
    const int INIT_NUM_OUTPUT_DUMPS = 5;
    const int INIT_NUM_CYCLES_PER_DUMP = 200;

    const double INIT_P_RATIO = 0.5;
    const double INIT_D_RATIO = 0.5;

    const bool verboseOutput = false;

    /**
     * \brief Simple representation of the mesh for this 1D example
     *
     * \details Mesh contains a set of elements and a set of faces between elements.
     *         It also contains three subsets: a single inflow; a single outflow element; all internal 'tube' elements. (added 5/2015)
     *         The mesh contains the relations from faces to elements and from tube elements to faces.
     *
     * \note (5/2015) We are currently missing an implicit constant grid relation.
     *
     * \note For current implementation with explicit static (constant) relations.
     * \note We are missing a nice way to set the relation elements.
     *       It should not have to be done explicitly in each user code -- especially in common use cases
     *       Idea: We could have a relationInverter function that takes a relation from sets A to B
     *             and generates a relation from set B to set A with all the arrows reversed.
     */
    class ShockTubeMesh
    {
    public:
        // types for sets
        typedef asctoolkit::meshapi::RangeSet ElemSet;
        typedef asctoolkit::meshapi::RangeSet FaceSet;

        // types for relations
        typedef asctoolkit::meshapi::StaticConstantRelation ElemToFaceRelation;
        typedef asctoolkit::meshapi::StaticConstantRelation FaceToElemRelation;

        // other types
        typedef asctoolkit::meshapi::Set::IndexType IndexType;
        typedef asctoolkit::meshapi::Set::PositionType PositionType;

    public:
        ElemSet elems;          // The entire set of elements
        ElemSet tubeElems;      // Subset of internal elements
        ElemSet inFlowElems;    // Subset of inflow elements (not used in this example)
        ElemSet outFlowElems;   // Subset of outflow elements (not used in this example)

        FaceSet faces;          // Faces between adjacent pairs of elements

        FaceToElemRelation relationFaceElem;    // co-boundary relation of faces to their elements
        ElemToFaceRelation relationTubeFace;    // boundary relation of internal 'tube' elements
    };


    // Define the explicit instances of our local (key/value) datastore for int and double
    // TODO: Might need an additional uint version for mesh data
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
   /* Get physics info */
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
      int numOutputDumps;
      int numCyclesPerDump ;

      std::cout << "How many dumps would you like? ";
      numOutputDumps = INIT_NUM_OUTPUT_DUMPS;
      std::cout << numOutputDumps <<std::endl;

      std::cout << "How many cycles between dumps would you like? ";
      numCyclesPerDump = INIT_NUM_CYCLES_PER_DUMP;
      std::cout << numCyclesPerDump <<std::endl;

      int numTotalCycles = numOutputDumps*numCyclesPerDump;
      std::cout << "\nSimulation will run for " << numTotalCycles <<" cycles."<<std::endl;

      intsRegistry.addScalar("numOutputDumps", numOutputDumps) ;
      intsRegistry.addScalar("numCyclesPerDump", numCyclesPerDump) ;
      intsRegistry.addScalar("numTotalCycles", numTotalCycles) ;
   }

   return ;
}

/**
 * \brief Build an empty mesh for the shock tube
 * \details Shocktube mesh layout
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

   // define the subsets
   ShockTubeMesh::PositionType numElems = mesh->elems.size();
   mesh->inFlowElems    = ShockTubeMesh::ElemSet( 0,1, &(mesh->elems) );
   mesh->tubeElems      = ShockTubeMesh::ElemSet( 1,numElems-1, &(mesh->elems) );
   mesh->outFlowElems   = ShockTubeMesh::ElemSet( numElems-1,numElems, &(mesh->elems) );

   // ------------ Set up relations

   // TODO: Need to define DynamicConstantRelation -- which will allow modifying the relation elements
   // TODO: Need to define ImplicitConstantRelation -- since the relations are actually implicit
   //       -- no storage should be necessary for regular grid neighbors
   // For now, we will have to do this explicitly...
   const ShockTubeMesh::PositionType STRIDE = 2;
   typedef std::vector<ShockTubeMesh::PositionType> IndexVec;

   /// Setup the FaceToElem relation
   IndexVec relVec( STRIDE * mesh->faces.size());
   IndexVec::iterator relIt = relVec.begin();
   for(ShockTubeMesh::IndexType idx =0; idx < static_cast<ShockTubeMesh::IndexType>(mesh->faces.size()); ++idx)
   {
       *relIt++ = idx;
       *relIt++ = idx+1;
   }

   mesh->relationFaceElem = ShockTubeMesh::FaceToElemRelation(&mesh->faces, &mesh->elems);
   mesh->relationFaceElem.bindRelationData(relVec, STRIDE);
   ATK_ASSERT(mesh->relationFaceElem.isValid( verboseOutput ));


   /// Setup the TubeElementToFace relation: A relation from the tubes subset of the elements to their incident faces
   //  For convenience, we can reuse the relVec container
   ShockTubeMesh::PositionType numTubeElems = mesh->tubeElems.size();
   relVec.clear();
   relVec.resize( STRIDE * numTubeElems);
   relIt = relVec.begin();
   for(ShockTubeMesh::IndexType idx =0; idx < static_cast<ShockTubeMesh::IndexType>(numTubeElems); ++idx)
   {
       *relIt++ = mesh->tubeElems[idx]-1;
       *relIt++ = mesh->tubeElems[idx];
   }

   mesh->relationTubeFace = ShockTubeMesh::ElemToFaceRelation(&mesh->tubeElems, &mesh->faces);
   mesh->relationTubeFace.bindRelationData(relVec, STRIDE);
   ATK_ASSERT(mesh->relationTubeFace.isValid( verboseOutput ));

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
   realsRegistry.addField("F0", &mesh.faces );      // mv
   realsRegistry.addField("F1", &mesh.faces );      // mv^2+P
   realsRegistry.addField("F2", &mesh.faces );      // v(E+P)

   // Fill left half with high pressure, right half with low pressure
   IndexType endTube    = mesh.elems.size();
   IndexType midTube    = endTube / 2;

   // Non-dimensionalized reference values
   double massInitial       = 1.0 ;
   double momentumInitial   = 0.0 ;
   double pressureInitial   = gammaaInverse ;
   double energyInitial     = pressureInitial/(gammaa-1.0) ;

   // Initialize zonal quantities
   RangeSet lowerTube(0, midTube);
   for (RangeSet::iterator elemIt=lowerTube.begin(); elemIt < lowerTube.end(); ++elemIt)
   {
      mass[*elemIt]     = massInitial ;
      momentum[*elemIt] = momentumInitial ;
      pressure[*elemIt] = pressureInitial ;
      energy[*elemIt]   = energyInitial ;
   }

   // adjust parameters for low pressure portion of tube
   massInitial      *= realsRegistry.getScalar("densityRatio") ;
   pressureInitial  *= realsRegistry.getScalar("pressureRatio") ;
   energyInitial     = pressureInitial/(gammaa - 1.0) ;

   RangeSet upperTube(midTube, mesh.elems.size());
   for (RangeSet::iterator elemIt=upperTube.begin(); elemIt < upperTube.end(); ++elemIt)
   {
       mass[*elemIt]     = massInitial ;
       momentum[*elemIt] = momentumInitial ;
       pressure[*elemIt] = pressureInitial ;
       energy[*elemIt]   = energyInitial ;
   }

   // Create needed time info
   realsRegistry.addScalar("time", 0.0) ;
   intsRegistry.addScalar("cycle", 0) ;

   double dx = 1.0 / static_cast<double>(endTube);
   realsRegistry.addScalar("dx", dx);
   realsRegistry.addScalar("dt", 0.4*dx) ;

}

/**
 * \brief Compute F quantities at faces.
 *
 * \details
 *
 * \verbatim
 *  @F   @F0   @F1   @F2
 *  -- = --- + --- + ---  
 *  @x   @x    @x    @x
 *  \endverbatim
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
   const RealField & mass     = realsRegistry.getField("mass") ;
   const RealField & momentum = realsRegistry.getField("momentum") ;
   const RealField & energy   = realsRegistry.getField("energy") ;

   // Update face data using element data using the face->elem relation
   ShockTubeMesh::PositionType numFaceElems = static_cast<ShockTubeMesh::PositionType>(mesh.faces.size() );
   for (ShockTubeMesh::PositionType fIdx = 0; fIdx < numFaceElems; ++fIdx)
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
   const RealField & F0 = realsRegistry.getField("F0") ;
   const RealField & F1 = realsRegistry.getField("F1") ;
   const RealField & F2 = realsRegistry.getField("F2") ;

   double   dx   = realsRegistry.getScalar("dx") ;
   double   dt   = realsRegistry.getScalar("dt") ;
   double & time = realsRegistry.getScalar("time") ;

   /// Update the element fields based on the face data using the elem->face relation
   ShockTubeMesh::PositionType numTubeElems = static_cast<ShockTubeMesh::PositionType>(mesh.tubeElems.size() );
   for (ShockTubeMesh::PositionType tPos=0; tPos < numTubeElems; ++tPos)
   {
      // Relation is over tube elements.
      ShockTubeMesh::IndexType elemIdx = mesh.tubeElems[tPos];

      // Each element inside the tube has an upwind and downwind face
      ShockTubeMesh::IndexType upWind   = mesh.relationTubeFace[tPos][UPWIND] ;      // upwind face
      ShockTubeMesh::IndexType downWind = mesh.relationTubeFace[tPos][DOWNWIND] ;    // downwind face

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

    static const ShockTubeMesh::PositionType MAX_ELEM_DUMP = 10;
    const int maxDump = std::min(mesh.elems.size(), MAX_ELEM_DUMP);
    const int rmaxDump = std::min(MAX_ELEM_DUMP, mesh.elems.size() - maxDump);

    // TODO: The following is currently grabbing the raw data from the Map and spitting out at most MAX_ELEM_DUMP elements
    // I would like to create a subset with a stride to only print every n_th element
    // Alternatively -- it can use an indirection map to grab the values, and write out to an sstream
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

   // We should be able to parallelize pretty easily by
   // adding an MPI_Init() here, modifying the setup slightly,
   // adding a communication subroutine in the main loop,
   // and calling MPI_Finalize() at the end of main()

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
