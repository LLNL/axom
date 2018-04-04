/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include "quest/quest.hpp"

// Quest includes
#include "axom/Types.hpp"
#include "axom/Macros.hpp"

#include "slic/slic.hpp"

#ifdef AXOM_USE_MPI
  #ifdef AXOM_USE_LUMBERJACK
    #include "slic/LumberjackStream.hpp"
  #else
    #include "slic/SynchronizedStream.hpp"
  #endif
#else
  #include "slic/GenericOutputStream.hpp"
#endif

#include "primal/BoundingBox.hpp"

#include "quest/InOutOctree.hpp"
#include "quest/SignedDistance.hpp"

#include "quest/STLReader.hpp"
#ifdef AXOM_USE_MPI
  #include "quest/PSTLReader.hpp"
#endif

namespace axom
{
namespace quest
{

// NOTE: supporting just one region/surface for now
// Note: Define everything in a local namespace
namespace
{

namespace slic = axom::slic;

typedef axom::mint::UnstructuredMesh< MINT_TRIANGLE > TriangleMesh;
enum QueryMode { QUERY_MODE_NONE,
                 QUERY_MODE_CONTAINMENT,
                 QUERY_MODE_SIGNED_DISTANCE };

/*!
 * \struct
 * \brief A simple struct to encapsulate knowledge about which
 *  acceleration structure we are using -- the SignedDistance or the InOutOctree
 */
template<int DIM>
struct QuestAccelerator
{
  typedef BoundingBox< double, DIM> GeometricBoundingBox;
  typedef Point< double, DIM> SpacePt;
  typedef Vector< double, DIM> SpaceVec;

  /** \brief Default constructor */
  QuestAccelerator()
    : m_surface_mesh(AXOM_NULLPTR),
    m_region(AXOM_NULLPTR),
    m_containmentTree(AXOM_NULLPTR),
    m_queryMode(QUERY_MODE_NONE),
    m_originalLoggerName(""),
    m_shouldFinalizeSlic(false),
    m_shouldDeleteMesh(true)
  {}

  /*!
   * \brief Sets the internal mesh pointer and computes some surface
   *  properties (bounding box and center of mass)
   */
  void setMesh( axom::mint::Mesh*& surface_mesh, bool deleteMesh )
  {
    SLIC_ASSERT( surface_mesh != AXOM_NULLPTR);

    m_surface_mesh = surface_mesh;
    m_shouldDeleteMesh = deleteMesh;

    // Compute the mesh's bounding box and center of mass
    m_meshBoundingBox.clear();
    m_meshCenterOfMass = SpacePt::zero();

    SpacePt pt;
    int numMeshNodes = m_surface_mesh->getNumberOfNodes();
    for ( int i=0 ; i < numMeshNodes ; ++i )
    {
      m_surface_mesh->getMeshNode( i, pt.data() );

      m_meshBoundingBox.addPoint( pt );
      m_meshCenterOfMass.array() += pt.array();
    }
    m_meshCenterOfMass.array() /= numMeshNodes;

    SLIC_ASSERT( m_meshBoundingBox.isValid() );
  }

  /*!
   * \brief Initializes the containment tree mode
   * \param surface_mesh The surface mesh
   * \pre Assumes that we are not yet initialized
   */
  void initializeContainmentTree( axom::mint::Mesh*& surface_mesh,
                                  bool deleteMesh )
  {
    SLIC_ASSERT( m_queryMode == QUERY_MODE_NONE);

    setMesh(surface_mesh, deleteMesh);
    m_containmentTree =
      new InOutOctree<DIM>( m_meshBoundingBox, m_surface_mesh );
    m_containmentTree->generateIndex();
    surface_mesh = m_surface_mesh;
    m_queryMode = QUERY_MODE_CONTAINMENT;
  }

  /*!
   * \brief Initializes the signed distance mode
   * \param surface_mesh The surface mesh
   * \pre Assumes that we are not yet initialized
   */
  void initializeSignedDistance( axom::mint::Mesh*& surface_mesh,
                                 int maxElements,
                                 int maxLevels,
                                 bool deleteMesh )
  {
    SLIC_ASSERT( m_queryMode == QUERY_MODE_NONE);

    setMesh(surface_mesh, deleteMesh);
    m_region =
      new SignedDistance<DIM>( m_surface_mesh, maxElements, maxLevels );
    surface_mesh = m_surface_mesh;
    m_queryMode = QUERY_MODE_SIGNED_DISTANCE;
  }

  /*!
   * \brief Deallocates all memory and sets the state to uninitialized
   */
  void finalize()
  {
    if ( m_region != AXOM_NULLPTR )
    {
      delete m_region;
      m_region = AXOM_NULLPTR;
    }

    if( m_containmentTree != AXOM_NULLPTR )
    {
      delete m_containmentTree;
      m_containmentTree = AXOM_NULLPTR;
    }
    m_queryMode = QUERY_MODE_NONE;

    if ( m_shouldDeleteMesh && m_surface_mesh != AXOM_NULLPTR )
    {

      delete m_surface_mesh;
      m_surface_mesh = AXOM_NULLPTR;
    }

    m_meshBoundingBox.clear();
    m_meshCenterOfMass = SpacePt::zero();
  }


  /*!
   * \brief Performs the distance query with the 3D point (x, y, z)
   * \param x The x-coordinate of the point
   * \param y The y-coordinate of the point
   * \param z The z-coordinate of the point
   * \return The signed distance from the point to the closest point on the
   * surface
   * \note Positive distances are outside the surface, negative distances are
   * inside.
   */
  double distance(double x, double y, double z)
  {
    SLIC_ASSERT_MSG(
      supportsDistanceQuery(),
      "Distance queries only supported when Quest is initialized with "
      << " requiresDistance = true." );

    SpacePt pt = SpacePt::make_point(x,y,z);
    return m_region->computeDistance( pt );
  }

  /*!
   * \brief Performs the containment query with the 3D point (x, y, z)
   * \param x The x-coordinate of the point
   * \param y The y-coordinate of the point
   * \param z The z-coordinate of the point
   * \return 1 if the point is inside the surface, 0 otherwise
   */
  int inside(double x, double y, double z)
  {
    SLIC_ASSERT( supportsContainmentQuery() );

    // TOOD: assume 3-D for now
    SpacePt pt = SpacePt::make_point(x,y,z);

    int sign = -1;
    switch( m_queryMode )
    {
    case QUERY_MODE_CONTAINMENT:
      sign = m_containmentTree->within(pt) ? 1 : 0;
      break;
    case QUERY_MODE_SIGNED_DISTANCE:
    {
      const quest::SignedDistance<3>::BVHTreeType* tree =
        m_region->getBVHTree();
      SLIC_ASSERT( tree != AXOM_NULLPTR );

      if ( !tree->contains( pt ) )
      {
        sign = 0;
      }
      else
      {
        sign = ( m_region->computeDistance( pt ) < 0.0f ) ? 1 : 0;
      }
    }
    break;
    case QUERY_MODE_NONE:
      break;
    }

    return( sign );
  }

  /*!
   * \brief Returns a const reference to the bounding box of the mesh
   */
  const GeometricBoundingBox& meshBoundingBox() const
  {
    return m_meshBoundingBox;
  }

  /*!
   * \brief Returns a const reference to the center of mass of the mesh
   */
  const SpacePt& meshCenterOfMass() const
  {
    return m_meshCenterOfMass;
  }

      #ifdef AXOM_DEBUG
  /*!
   * \brief Utility function to determine if we are in a mode that supports
   * distance queries
   */
  bool supportsDistanceQuery()
  {
    bool isValid = true;

    switch(m_queryMode)
    {
    case QUERY_MODE_CONTAINMENT:
      isValid = false;
      break;
    case QUERY_MODE_SIGNED_DISTANCE:
      if( m_region == AXOM_NULLPTR)
        isValid = false;
      break;
    case QUERY_MODE_NONE:
      isValid = false;
      break;
    }

    return isValid;
  }

  /*!
   * \brief Utility function to determine if we are in a mode that supports
   * containment queries
   */
  bool supportsContainmentQuery()
  {
    bool isValid = true;

    switch(m_queryMode)
    {
    case QUERY_MODE_CONTAINMENT:
      if( m_containmentTree == AXOM_NULLPTR)
        isValid = false;
      break;
    case QUERY_MODE_SIGNED_DISTANCE:
      if( m_region == AXOM_NULLPTR)
        isValid = false;
      break;
    case QUERY_MODE_NONE:
      isValid = false;
      break;
    }

    return isValid;
  }

  /*!
   * \brief Utility function to determine if an acceleration structure has been
   * initialized
   */
  bool isInitialized()
  {
    return m_queryMode == QUERY_MODE_CONTAINMENT ||
           m_queryMode == QUERY_MODE_SIGNED_DISTANCE;
  }
  #endif

  /*!
   * \brief Sets up the formatted Slic logger for quest
   */
#ifdef AXOM_USE_MPI
  void setupQuestLogger( MPI_Comm comm)
#else
  void setupQuestLogger()
#endif
  {
    // Ensure slic has been initialized
    if( !slic::isInitialized() )
    {
      slic::initialize();
      m_shouldFinalizeSlic = true;            // mark that we need to finalize
                                              // slic
    }

    const std::string questLoggerName = "quest_logger";
    m_originalLoggerName = slic::getActiveLoggerName();

    slic::flushStreams();
    if( !slic::activateLogger(questLoggerName) )
    {
      slic::LogStream* ls;

              #ifdef AXOM_USE_MPI
      std::string fmt = "[<RANK>][Quest <LEVEL>]: <MESSAGE>\n";
                #ifdef AXOM_USE_LUMBERJACK
      const int RLIMIT = 8;
      ls = new slic::LumberjackStream( &std::cout, comm, RLIMIT, fmt);
                #else
      ls = new slic::SynchronizedStream( &std::cout, comm, fmt );
                #endif
              #else
      std::string fmt = "[Quest <LEVEL>]: <MESSAGE>\n";
      ls = new slic::GenericOutputStream(&std::cout, fmt);
              #endif
      slic::createLogger(questLoggerName);
      slic::activateLogger(questLoggerName);
      slic::setLoggingMsgLevel( slic::message::Info );
      slic::addStreamToAllMsgLevels(ls);
    }
  }

  /*!
   * \brief Deactivates the quest logger
   *
   * If there was a previous logger, it is restored.
   * If slic was initialized by a call to setupQuestLogger(), slic is finalized.
   */
  void teardownQuestLogger()
  {
    if( !slic::isInitialized())
    {
      m_shouldFinalizeSlic = false;
      return;
    }

    slic::flushStreams();

    if(m_originalLoggerName != "")
    {
      // Revert to original Slic logger
      slic::activateLogger(m_originalLoggerName);
      m_originalLoggerName = "";
    }

    if( m_shouldFinalizeSlic )
    {
      slic::finalize();
      m_shouldFinalizeSlic = false;
    }

  }


private:
  axom::mint::Mesh* m_surface_mesh;
  SignedDistance< DIM >* m_region;
  InOutOctree< DIM >* m_containmentTree;
  QueryMode m_queryMode;

  SpacePt m_meshCenterOfMass;
  GeometricBoundingBox m_meshBoundingBox;

  std::string m_originalLoggerName;
  bool m_shouldFinalizeSlic;
  bool m_shouldDeleteMesh;
};

/*!
 * \brief A static instance of the acceleration structure in 3D
 * \note In this initial release, we assume a single static accelerator.
 *  Eventually, we will expand on this to support multiple structures in 2D
 *  and 3D. We will probably use Sidre to hold pointers to these structures.
 */
static QuestAccelerator<3> accelerator3D;
}


//------------------------------------------------------------------------------
#ifdef AXOM_USE_MPI
void initialize( MPI_Comm comm, const std::string& fileName,
                 bool requiresDistance, int ndims, int maxElements,
                 int maxLevels )
{
  SLIC_ASSERT( comm != MPI_COMM_NULL );

  // Read in the mesh
  quest::PSTLReader* reader = new quest::PSTLReader( comm );
  reader->setFileName( fileName );
  reader->read();

  axom::mint::Mesh* surface_mesh = new TriangleMesh( 3 );
  SLIC_ASSERT( surface_mesh != AXOM_NULLPTR );

  reader->getMesh( static_cast< TriangleMesh* >( surface_mesh ) );
  delete reader;

  initialize(comm, surface_mesh, requiresDistance, ndims, maxElements,
             maxLevels, true);
}

void initialize( MPI_Comm comm, mint::Mesh*& input_mesh,
                 bool requiresDistance, int ndims, int maxElements,
                 int maxLevels, bool deleteMesh )
{
  SLIC_ASSERT( !accelerator3D.isInitialized() );
  SLIC_ASSERT( comm != MPI_COMM_NULL );

  AXOM_DEBUG_VAR(ndims);
  SLIC_ASSERT( ndims==2 || ndims==3 );

  // In the future, we will also support 2D, but we currently only support 3D
  SLIC_ASSERT_MSG(ndims==3,
                  "Quest currently only supports 3D (not 2D) triangle meshes.");
  SLIC_ASSERT_MSG(input_mesh->getMeshType() == MINT_UNSTRUCTURED_TRIANGLE_MESH,
                  "Quest currently only supports 3D triangle meshes "
                  "(not any other kind of cell).");

  accelerator3D.setupQuestLogger(comm);

  SLIC_ASSERT( input_mesh != AXOM_NULLPTR );

  // Initialize the appropriate acceleration structure
  if(requiresDistance)
  {
    accelerator3D.initializeSignedDistance(input_mesh, maxElements,
                                           maxLevels, deleteMesh);
  }
  else
  {
    accelerator3D.initializeContainmentTree(input_mesh, deleteMesh);
  }

  accelerator3D.teardownQuestLogger();
}
#else
//------------------------------------------------------------------------------
void initialize( const std::string& fileName,
                 bool requiresDistance, int ndims, int maxElements,
                 int maxLevels )
{
  // Read in the mesh
  quest::STLReader* reader = new quest::STLReader();
  reader->setFileName( fileName );
  reader->read();

  axom::mint::Mesh* surface_mesh = new TriangleMesh( 3 );
  SLIC_ASSERT( surface_mesh != AXOM_NULLPTR );

  reader->getMesh( static_cast< TriangleMesh* >( surface_mesh ) );
  delete reader;

  initialize(surface_mesh, requiresDistance, ndims, maxElements,
             maxLevels, true);
}

void initialize( mint::Mesh*& input_mesh,
                 bool requiresDistance, int ndims, int maxElements,
                 int maxLevels, bool deleteMesh )
{
  SLIC_ASSERT( !accelerator3D.isInitialized() );

  AXOM_DEBUG_VAR(ndims);
  SLIC_ASSERT( ndims==2 || ndims==3 );

  // In the future, we will also support 2D, but we currently only support 3D
  SLIC_ASSERT_MSG(ndims==3,
                  "Quest currently only supports 3D (not 2D) triangle meshes.");
  SLIC_ASSERT_MSG(input_mesh->getMeshType() == MINT_UNSTRUCTURED_TRIANGLE_MESH,
                  "Quest currently only supports 3D triangle meshes "
                  "(not any other kind of cell).");

  accelerator3D.setupQuestLogger();

  SLIC_ASSERT( input_mesh != AXOM_NULLPTR );

  // Initialize the appropriate acceleration structure
  if(requiresDistance)
  {
    accelerator3D.initializeSignedDistance(input_mesh, maxElements,
                                           maxLevels, deleteMesh );
  }
  else
  {
    accelerator3D.initializeContainmentTree(input_mesh, deleteMesh);
  }

  accelerator3D.teardownQuestLogger();
}
#endif

//------------------------------------------------------------------------------
double distance( double x, double y, double z )
{
  // TODO: assume 3-D for now
  return accelerator3D.distance(x,y,z);
}

//------------------------------------------------------------------------------
void distance( const double* xyz, double* dist, int npoints )
{
  SLIC_ASSERT( xyz != AXOM_NULLPTR );
  SLIC_ASSERT( dist != AXOM_NULLPTR );

#ifdef AXOM_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for ( int i=0 ; i < npoints ; ++i )
  {
    // TODO: assume 3-D for now
    dist[ i ] = accelerator3D.distance(xyz[i*3], xyz[i*3+1], xyz[i*3+2]);
  }

}

//------------------------------------------------------------------------------
int inside( double x, double y, double z )
{
  // TODO: assume 3-D for now
  return accelerator3D.inside(x,y,z);
}




//------------------------------------------------------------------------------
void mesh_min_bounds(double* coords)
{
  typedef QuestAccelerator<3>::SpacePt SpacePt;
  SLIC_ASSERT(coords != AXOM_NULLPTR);

  const SpacePt& bbMin = accelerator3D.meshBoundingBox().getMin();
  bbMin.array().to_array(coords);
}

//------------------------------------------------------------------------------
void mesh_max_bounds(double* coords)
{
  typedef QuestAccelerator<3>::SpacePt SpacePt;
  SLIC_ASSERT(coords != AXOM_NULLPTR);

  const SpacePt& bbMax = accelerator3D.meshBoundingBox().getMax();
  bbMax.array().to_array(coords);
}



//------------------------------------------------------------------------------
void mesh_center_of_mass(double* coords)
{
  typedef QuestAccelerator<3>::SpacePt SpacePt;
  SLIC_ASSERT(coords != AXOM_NULLPTR);

  const SpacePt& cMass = accelerator3D.meshCenterOfMass();
  cMass.array().to_array(coords);
}

//------------------------------------------------------------------------------
void inside( const double* xyz, int* in, int npoints )
{
  SLIC_ASSERT( xyz != AXOM_NULLPTR );
  SLIC_ASSERT( in != AXOM_NULLPTR );

#ifdef AXOM_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for ( int i=0 ; i < npoints ; ++i )
  {
    // TODO: assume 3-D for now
    in[ i ] = accelerator3D.inside( xyz[i*3], xyz[i*3+1], xyz[i*3+2] );
  }
}

//------------------------------------------------------------------------------
void finalize()
{
  accelerator3D.finalize();
  accelerator3D.teardownQuestLogger();
}

} // end namespace quest
} // end namespace axom
