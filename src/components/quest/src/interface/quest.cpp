/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

#include "quest/quest.hpp"

// Quest includes
#include "common/CommonTypes.hpp"
#include "quest/STLReader.hpp"
#include "quest/SignedDistance.hpp"

#include "quest/BoundingBox.hpp"
#include "quest/InOutOctree.hpp"

#include "slic/slic.hpp"

#ifdef USE_MPI
#include "quest/PSTLReader.hpp"
#endif

namespace quest
{
    // NOTE: supporting just one region/surface for now
    // Note: Define everything in a local namespace
  namespace
  {
    typedef mint::UnstructuredMesh< mint::LINEAR_TRIANGLE > TriangleMesh;
    enum QueryMode { QUERY_MODE_NONE, QUERY_MODE_CONTAINMENT, QUERY_MODE_SIGNED_DISTANCE };

    /**
     * \struct
     * \brief A simple struct to encapsulate knowledge about which acceleration structure we are using
     *        -- the SignedDistance or the InOutOctree
     */
    template<int DIM>
    struct QuestAccelerator
    {
        typedef BoundingBox< double, DIM> GeometricBoundingBox;
        typedef Point< double, DIM> SpacePt;
        typedef Vector< double, DIM> SpaceVec;

        /** \brief Default constructor */
        QuestAccelerator()
            : m_surface_mesh(ATK_NULLPTR)
            , m_region(ATK_NULLPTR)
            , m_containmentTree(ATK_NULLPTR)
            , m_queryMode(QUERY_MODE_NONE) {}

        /**
         * \brief Sets the internal mesh pointer and computes some surface properties (bounding box and center of mass)
         */
        void setMesh(mint::Mesh* surface_mesh)
        {
            SLIC_ASSERT( surface_mesh != ATK_NULLPTR);

            m_surface_mesh = surface_mesh;

            // Computer the mesh's bounding box and center of mass
            m_meshBoundingBox.clear();
            m_meshCenterOfMass = SpacePt::zero();

            SpacePt pt;
            int numMeshNodes = m_surface_mesh->getMeshNumberOfNodes();
            for ( int i=0; i < numMeshNodes; ++i )
            {
                m_surface_mesh->getMeshNode( i, pt.data() );

                m_meshBoundingBox.addPoint( pt );
                m_meshCenterOfMass.array() += pt.array();
            }
            m_meshCenterOfMass.array() /= numMeshNodes;

            SLIC_ASSERT( m_meshBoundingBox.isValid() );
        }

        /**
         * \brief Initializes the containment tree mode
         * \param surface_mesh The surface mesh
         * \pre Assumes that we are not yet initialized
         */
        void initializeContainmentTree(mint::Mesh* surface_mesh)
        {
            SLIC_ASSERT( m_queryMode == QUERY_MODE_NONE);

            setMesh(surface_mesh);
            m_containmentTree = new InOutOctree<DIM>( m_meshBoundingBox, m_surface_mesh );
            m_containmentTree->generateIndex();
            m_queryMode = QUERY_MODE_CONTAINMENT;
        }

        /**
         * \brief Initializes the signed distance mode
         * \param surface_mesh The surface mesh
         * \pre Assumes that we are not yet initialized
         */
        void initializeSignedDistance(mint::Mesh* surface_mesh, int maxElements, int maxLevels)
        {
            SLIC_ASSERT( m_queryMode == QUERY_MODE_NONE);

            setMesh(surface_mesh);
            m_region = new SignedDistance<DIM>( m_surface_mesh, maxElements, maxLevels );
            m_queryMode = QUERY_MODE_SIGNED_DISTANCE;
        }

        /**
         * \brief Deallocates all memory and sets the state to uninitialized
         */
        void finalize()
        {
            if ( m_region != ATK_NULLPTR ) {
               delete m_region;
               m_region = ATK_NULLPTR;
            }

            if( m_containmentTree != ATK_NULLPTR )
            {
                delete m_containmentTree;
                m_containmentTree = ATK_NULLPTR;
            }
            m_queryMode = QUERY_MODE_NONE;

            if ( m_surface_mesh != ATK_NULLPTR ) {

               delete m_surface_mesh;
               m_surface_mesh = ATK_NULLPTR;
            }

            m_meshBoundingBox.clear();
            m_meshCenterOfMass = SpacePt::zero();
        }


        /**
         * \brief Performs the distance query with the 3D point (x, y, z)
         * \param x The x-coordinate of the point
         * \param y The y-coordinate of the point
         * \param z The z-coordinate of the point
         * \return The signed distance from the point to the closest point on the surface
         * \note Positive distances are outside the surface, negative distances are inside.
         */
        double distance(double x, double y, double z)
        {
            SLIC_ASSERT_MSG( supportsDistanceQuery()
                           , "Distance queries only supported when Quest is initialized with "
                             << " requiresDistance = true." );

            SpacePt pt = SpacePt::make_point(x,y,z);
            return m_region->computeDistance( pt );
        }

        /**
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
                const quest::SignedDistance<3>::BVHTreeType* tree = m_region->getBVHTree();
                SLIC_ASSERT( tree != ATK_NULLPTR );

                if ( !tree->contains( pt ) ) {
                  sign = 0;
                } else {
                  sign = ( m_region->computeDistance( pt ) < 0.0f )? 1 : 0;
                }
              }
              break;
            case QUERY_MODE_NONE:
              break;
            }

            return( sign );
        }

        /**
         * \brief Returns a const reference to the bounding box of the mesh
         */
        const GeometricBoundingBox& meshBoundingBox() const
        {
            return m_meshBoundingBox;
        }

        /**
         * \brief Returns a const reference to the center of mass of the mesh
         */
        const SpacePt& meshCenterOfMass() const
        {
            return m_meshCenterOfMass;
        }

      #ifdef ATK_DEBUG
        /**
         * \brief Utility function to determine if we are in a mode that supports distance queries
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
                if( m_region == ATK_NULLPTR)
                    isValid = false;
                break;
            case QUERY_MODE_NONE:
                isValid = false;
                break;
            }

            return isValid;
        }

        /**
         * \brief Utility function to determine if we are in a mode that supports containment queries
         */
        bool supportsContainmentQuery()
        {
            bool isValid = true;

            switch(m_queryMode)
            {
            case QUERY_MODE_CONTAINMENT:
                if( m_containmentTree == ATK_NULLPTR)
                    isValid = false;
                break;
            case QUERY_MODE_SIGNED_DISTANCE:
                if( m_region == ATK_NULLPTR)
                    isValid = false;
                break;
            case QUERY_MODE_NONE:
                isValid = false;
                break;
            }

            return isValid;
        }

        /**
         * \brief Utility function to determine if an acceleration structure has been initialized
         */
        bool isInitialized()
        {
            return m_queryMode == QUERY_MODE_CONTAINMENT || m_queryMode == QUERY_MODE_SIGNED_DISTANCE;
        }
  #endif

    private:
        mint::Mesh* m_surface_mesh;
        SignedDistance< DIM >* m_region;
        InOutOctree< DIM >* m_containmentTree;
        QueryMode m_queryMode;

        SpacePt              m_meshCenterOfMass;
        GeometricBoundingBox m_meshBoundingBox;
    };

    /**
     * \brief A static instance of the acceleration structure in 3D
     * \note In this initial release, we assume a single static accelerator.
     *       Eventually, we will expand on this to support multiple structures in 2D and 3D.
     *       We will probably use Sidre to hold pointers to these structures.
     */
    static QuestAccelerator<3> accelerator3D;
}


//------------------------------------------------------------------------------
#ifdef USE_MPI
void initialize( MPI_Comm comm, const std::string& fileName,
                 bool requiresDistance, int ndims, int maxElements, int maxLevels )
{
  SLIC_ASSERT( ! accelerator3D.isInitialized() );
  SLIC_ASSERT( comm != MPI_COMM_NULL );

  ATK_DEBUG_VAR(ndims);
  SLIC_ASSERT( ndims==2 || ndims==3 );

  // In the future, we will also support 2D, but we currently only support 3D
  SLIC_ASSERT_MSG(ndims==3, "Quest currently only supports 3D triangle meshes.");

  // Read in the mesh
  quest::PSTLReader* reader = new quest::PSTLReader( comm );
  reader->setFileName( fileName );
  reader->read();

  mint::Mesh* surface_mesh = new TriangleMesh( 3 );
  SLIC_ASSERT( surface_mesh != ATK_NULLPTR );

  reader->getMesh( static_cast< TriangleMesh* >( surface_mesh ) );
  delete reader;

  // Initialize the appropriate acceleration structure
  if(requiresDistance)
  {
      accelerator3D.initializeSignedDistance(surface_mesh, maxElements, maxLevels);
  }
  else
  {
      accelerator3D.initializeContainmentTree(surface_mesh);
  }
}
#else
//------------------------------------------------------------------------------
void initialize( const std::string& fileName,
                 bool requireDistance, int ndims, int maxElements, int maxLevels )
{
  SLIC_ASSERT( ! accelerator3D.isInitialized() );

  int NDIMS = ndims;
  SLIC_ASSERT( NDIMS==2 || NDIMS==3 );

  // In the future, we will also support 2D, but we currently only support 3D
  SLIC_ASSERT_MSG(NDIMS==3, "Quest currently only supports 3D triangle meshes.");

  // Read in the mesh
  quest::STLReader* reader = new quest::STLReader();
  reader->setFileName( fileName );
  reader->read();

  mint::Mesh* surface_mesh = new TriangleMesh( 3 );
  SLIC_ASSERT( surface_mesh != ATK_NULLPTR );

  reader->getMesh( static_cast< TriangleMesh* >( surface_mesh ) );
  delete reader;

  // Initialize the appropriate acceleration structure
  if(requireDistance)
  {
    accelerator3D.initializeSignedDistance(surface_mesh, maxElements, maxLevels );
  }
  else
  {
    accelerator3D.initializeContainmentTree(surface_mesh);
  }
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
  SLIC_ASSERT( xyz != ATK_NULLPTR );
  SLIC_ASSERT( dist != ATK_NULLPTR );

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for ( int i=0; i < npoints; ++i ) {
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
double mesh_min_x(){ return accelerator3D.meshBoundingBox().getMin()[0]; }

//------------------------------------------------------------------------------
double mesh_max_x(){ return accelerator3D.meshBoundingBox().getMax()[0]; }

//------------------------------------------------------------------------------
double mesh_min_y(){ return accelerator3D.meshBoundingBox().getMin()[1]; }

//------------------------------------------------------------------------------
double mesh_max_y(){ return accelerator3D.meshBoundingBox().getMax()[1]; }

//------------------------------------------------------------------------------
double mesh_min_z(){ return accelerator3D.meshBoundingBox().getMin()[2]; }

//------------------------------------------------------------------------------
double mesh_max_z(){ return accelerator3D.meshBoundingBox().getMax()[2]; }


//------------------------------------------------------------------------------
double mesh_center_of_mass_x(){ return accelerator3D.meshCenterOfMass()[0]; }

//------------------------------------------------------------------------------
double mesh_center_of_mass_y(){ return accelerator3D.meshCenterOfMass()[1]; }

//------------------------------------------------------------------------------
double mesh_center_of_mass_z(){ return accelerator3D.meshCenterOfMass()[2]; }



//------------------------------------------------------------------------------
void inside( const double* xyz, int* in, int npoints )
{
  SLIC_ASSERT( xyz != ATK_NULLPTR );
  SLIC_ASSERT( in != ATK_NULLPTR );

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for ( int i=0; i < npoints; ++i ) {
     // TODO: assume 3-D for now
     in[ i ] = accelerator3D.inside( xyz[i*3], xyz[i*3+1], xyz[i*3+2] );
  }
}

//------------------------------------------------------------------------------
void finalize()
{
  accelerator3D.finalize();
}

} /* end namespace quest */

