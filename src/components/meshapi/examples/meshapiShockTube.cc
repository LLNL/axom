/**
 * \file meshapiShockTube.cc
 *
 * \brief 1D shock tube, split flux Euler equations
 *
 * \author J. Keasler (original)
 * \author K. Weiss (modified to use the asc toolkit mesh API)
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


#include "common/Utilities.hpp"
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
    const int INIT_NUM_ULTRA_DUMPS = 10;
    const int INIT_NUM_ULTRA_CYCLES_PER_DUMP = 100;

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
        typedef asctoolkit::meshapi::StaticConstantRelation ElemFaceRelation;
        typedef asctoolkit::meshapi::StaticConstantRelation FaceElemRelation;

        // other types
        typedef asctoolkit::meshapi::MeshIndexType IndexType;

    public:
        ElemSet elems;
        ElemSet faces;

        FaceElemRelation relationFaceElem;
        ElemFaceRelation relationElemFace;

    };

    /** \brief A helper class to print the name of a few types */
    template<typename T> struct TypeToString{};

    /** \brief A helper class to print the name of integers as 'int' */
    template<> struct TypeToString<int>{ static std::string to_string(){return "int";} };

    /** \brief A helper class to print the name of doubles as 'double' */
    template<> struct TypeToString<double>{ static std::string to_string(){return "double";} };

    /**
     * \brief Very simple container for fields of a given type DataType with minimal error checking.
     * \note
     *         We are using concrete instances for int and double in the code below.
     *         This should eventually be replaced with the sidre datastore.
     */
    template<typename DataType>
    class FieldRegistry
    {
    public:
        typedef std::string                     KeyType;
        typedef typename std::vector<DataType>  DataVec;
        typedef std::map<KeyType, DataVec>      DataVecMap;
        typedef std::map<KeyType, DataType>     DataAttrMap;

    public:
        DataVec&  addField(KeyType key, DataVec& vec) { return m_dataVecs[key] = vec; }
        DataType& addScalar(KeyType key, DataType val) { return m_dataScalars[key] = val; }

        DataVec& getField(KeyType key)
        {
            ATK_ASSERT_MSG( m_dataVecs.find(key) != m_dataVecs.end(), "Didn't find " << TypeToString<DataType>::to_string() << " field named " << key );
            return m_dataVecs[key];
        }
        DataVec const& getField(KeyType key) const
        {
            ATK_ASSERT_MSG( m_dataVecs.find(key) != m_dataVecs.end(), "Didn't find " << TypeToString<DataType>::to_string() << " field named " << key );
            return m_dataVecs[key];
        }

        DataType& getScalar(KeyType key)
        {
            ATK_ASSERT_MSG( m_dataScalars.find(key) != m_dataScalars.end(), "Didn't find " << TypeToString<DataType>::to_string() << " scalar named " << key );
            return m_dataScalars[key];
        }
        DataType const& getScalar(KeyType key) const
        {
            ATK_ASSERT_MSG( m_dataScalars.find(key) != m_dataScalars.end(), "Didn't find " << TypeToString<DataType>::to_string() << " scalar named " << key );
            return m_dataScalars[key];
        }

    private:
        DataVecMap  m_dataVecs;
        DataAttrMap m_dataScalars;
    };

    // Define the explicit instances for int and double
    FieldRegistry<int>    intsRegistry;
    FieldRegistry<double> realsRegistry;

    typedef FieldRegistry<double>::DataVec RealField;


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
   for(ShockTubeMesh::IndexType idx =0; idx < mesh->faces.size(); ++idx)
   {
       *relIt++ = idx;
       *relIt++ = idx+1;
   }

   mesh->relationFaceElem = ShockTubeMesh::FaceElemRelation(&mesh->faces, &mesh->elems);
   mesh->relationFaceElem.setRelation(relVec, STRIDE);
   ATK_ASSERT(mesh->relationFaceElem.isValid());


   /// Setup the elem -> face relation
   // Note: This relation needs to be on the TUBE subset of the elements (i.e. skip the first and last elt
   //       It is currently on all elements, which necessitates a workaround below.
   //       And the first and last elements of the relation are set to 0 which is wrong.
   //       Question -- how are we planning to handle indexes that are out or range (accidentally)?
   //                   how are we planning to handle indexes that are intentionally out of range
   //                   (e.g. to indicate a problem, or a missing element etc..)?
   unsigned int const elemSize = mesh->elems.size();
   relVec.clear();
   relVec.resize( STRIDE * elemSize);
   relIt = relVec.begin();
   for(ShockTubeMesh::IndexType idx =0; idx < elemSize; ++idx)
   {
       // HACK -- these indexes should be over the tube subset -- skipping the first and last elem
       *relIt++ = (idx == 0)? 0 : idx-1;
       *relIt++ = (idx == elemSize-1)? 0 : idx;
   }

   mesh->relationElemFace = ShockTubeMesh::ElemFaceRelation(&mesh->elems, &mesh->faces);
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
   RealField elemField(mesh.elems.size());
   RealField& mass      = realsRegistry.addField("mass", elemField );
   RealField& momentum  = realsRegistry.addField("momentum", elemField );
   RealField& energy    = realsRegistry.addField("energy", elemField );
   RealField& pressure  = realsRegistry.addField("pressure", elemField );

   // Create face centered fields
   // mv, mv^2+P, and v(E+P)
   RealField faceField(mesh.faces.size());
   realsRegistry.addField("F0", faceField );
   realsRegistry.addField("F1", faceField );
   realsRegistry.addField("F2", faceField );

   // Fill left half with high pressure, right half with low pressure
   IndexType startTube = 0 ;
   IndexType endTube = mesh.elems.size();
   IndexType midTube = endTube / 2 ;

   // Non-dimensionalized reference values
   double massInitial = 1.0 ;
   double momentumInitial = 0.0 ;
   double pressureInitial = gammaaInverse ;
   double energyInitial = pressureInitial/(gammaa-1.0) ;

   // Initialize zonal quantities
   for (IndexType idx=startTube; idx<midTube; ++idx)
   {
      mass[idx] = massInitial ;
      momentum[idx] = momentumInitial ;
      pressure[idx] = pressureInitial ;
      energy[idx] = energyInitial ;
   }

   // adjust parameters for low pressure portion of tube
   double dratio = realsRegistry.getScalar("densityRatio") ;
   double pratio = realsRegistry.getScalar("pressureRatio") ;

   massInitial *= dratio ;
   pressureInitial *= pratio ;
   energyInitial = pressureInitial/(gammaa - 1.0) ;

   for (IndexType idx=midTube; idx<endTube; ++idx)
   {
      mass[idx] = massInitial ;
      momentum[idx] = momentumInitial ;
      pressure[idx] = pressureInitial ;
      energy[idx] = energyInitial ;
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
   RealField& F0 = realsRegistry.getField("F0") ;
   RealField& F1 = realsRegistry.getField("F1") ;
   RealField& F2 = realsRegistry.getField("F2") ;

   // Element fields
   RealField const& mass =     realsRegistry.getField("mass") ;
   RealField const& momentum = realsRegistry.getField("momentum") ;
   RealField const& energy =   realsRegistry.getField("energy") ;

   // Update face data using element data using the face->elem relation
   for (ShockTubeMesh::IndexType fIdx=0; fIdx< mesh.faces.size() ; ++fIdx)
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
      ev = v*(gammaa - 1.0) ;

      F0[fIdx] += ev*massf ;
      F1[fIdx] += ev*momentumf ;
      F2[fIdx] += ev*(energyf - pressuref) ;

      contributor = ((v + c >= 0.0) ? upWind : downWind) ;
      massf     = mass[contributor] ;
      momentumf = momentum[contributor] ;
      energyf   = energy[contributor] ;
      pressuref = (gammaa - 1.0)*(energyf - 0.5*momentumf*momentumf/massf) ;
      ev = 0.5*(v + c) ;
      cLocal = sqrt(gammaa*pressuref/massf) ;

      F0[fIdx] += ev*massf ;
      F1[fIdx] += ev*(momentumf + massf*cLocal) ;
      F2[fIdx] += ev*(energyf + pressuref + momentumf*cLocal) ;

      contributor = ((v - c >= 0.0) ? upWind : downWind) ;
      massf     = mass[contributor] ;
      momentumf = momentum[contributor] ;
      energyf   = energy[contributor] ;
      pressuref = (gammaa - 1.0)*(energyf - 0.5*momentumf*momentumf/massf) ;
      ev = 0.5*(v - c) ;
      cLocal = sqrt(gammaa*pressuref/massf) ;

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
    RealField& mass = realsRegistry.getField("mass") ;
    RealField& momentum = realsRegistry.getField("momentum") ;
    RealField& energy = realsRegistry.getField("energy") ;
    RealField& pressure = realsRegistry.getField("pressure") ;

   // The element update is calculated as the flux between faces
   RealField const& F0 = realsRegistry.getField("F0") ;
   RealField const& F1 = realsRegistry.getField("F1") ;
   RealField const& F2 = realsRegistry.getField("F2") ;

   double dx = realsRegistry.getScalar("dx") ;
   double dt = realsRegistry.getScalar("dt") ;
   double &time = realsRegistry.getScalar("time") ;


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
            //   DumpUltra(problem) ;
            std::cout<< "\n\tStarting cycle " << currCycle
                 << " at time " << realsRegistry.getScalar("time");
        }


        ComputeFaceInfo(mesh) ;
        UpdateElemInfo(mesh) ;
    }

    std::cout<<"\ndone." << std::endl;

    //DumpUltra(problem) ; /* One last dump */

    return 0 ;
}
