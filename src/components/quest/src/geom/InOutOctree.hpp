/**
 * \file
 *
 * \brief Defines an InOutOctree for containment queries on a surface.
 */

#ifndef INOUT_OCTREE__HXX_
#define INOUT_OCTREE__HXX_

#include "common/config.hpp"       // defines AXOM_USE_CXX11
#include "common/AxomMacros.hpp"
#include "common/Timer.hpp"
#include "common/Utilities.hpp"

#include "primal/BoundingBox.hpp"
#include "primal/Point.hpp"
#include "primal/Triangle.hpp"
#include "primal/Vector.hpp"

#include "primal/intersection.hpp"
#include "primal/orientation.hpp"
#include "primal/squared_distance.hpp"


#include "slic/slic.hpp"

#include "slam/Map.hpp"
#include "slam/RangeSet.hpp"
#include "slam/StaticConstantRelation.hpp"
#include "slam/StaticVariableRelation.hpp"

#include "quest/SpatialOctree.hpp"

#include "mint/Mesh.hpp"
#include "mint/UnstructuredMesh.hpp"
#include "mint/FieldData.hpp"
#include "mint/FieldVariable.hpp"


#include <vector>   // For InOutLeafData triangle lists -- TODO replace with SLAM DynamicVariableRelation...
#include <iterator> // For back_inserter
#include <limits>   // numeric_limits traits
#include <sstream>


#define DEBUG_VERT_IDX  -2
#define DEBUG_TRI_IDX   -2

#define DEBUG_BLOCK_2  BlockIndex::invalid_index()
#define DEBUG_BLOCK_1  BlockIndex::invalid_index()

#ifndef DUMP_VTK_MESH
//    #define DUMP_VTK_MESH
#endif

#ifndef DUMP_OCTREE_INFO
//    #define DUMP_OCTREE_INFO 1
#endif

#ifndef DEBUG_OCTREE_ACTIVE
//    #define DEBUG_OCTREE_ACTIVE
#endif

#if defined(DEBUG_OCTREE_ACTIVE) and defined(AXOM_DEBUG)
  #define QUEST_OCTREE_DEBUG_LOG_IF( _cond, _msg)  if( _cond ) SLIC_DEBUG( _msg )
#else
  #define QUEST_OCTREE_DEBUG_LOG_IF( _cond, _msg)  ( (void)0 )
#endif


namespace axom {  
namespace quest {

    /**
     * \brief Compact BlockDataType for an InOutOctree
     *
     * Storage requirement is one integer per block, which can implicitly determine
     * the color of a block, and for gray block, the index of the associated triangles
     */
    class InOutBlockData
    {
        // Some internal constants for keeping tracking of the associated block
        // A block is a leaf block when its m_idx is not INTERNAL_BLOCK
        // Leaf blocks can be uncolored or colored (without additional data)
        //      or m_idx be the index of the data associated with a gray block
        enum { LEAF_BLOCK_UNCOLORED = -1
             , LEAF_BLOCK_WHITE     = -2
             , LEAF_BLOCK_BLACK     = -3
             , INTERNAL_BLOCK       = -4
             , NON_BLOCK            = -5
        };

    public:
        enum LeafColor { Undetermined = -2, White = -1, Gray = 0, Black = 1};

    public:
        /**
         * \brief Default constructor for an InOutBlockData
         *
         * \note Default constructed InOutBlockData instances are assumed to be leaf blocks
         */
        InOutBlockData()
            : m_idx(LEAF_BLOCK_UNCOLORED) {}

        /** \brief Constructor from a given index */
        explicit InOutBlockData(int dataIdx)
            : m_idx(dataIdx){}

        /** \brief Copy constructor for an InOutBlockData instance */
        InOutBlockData(const InOutBlockData& other)
            : m_idx(other.m_idx){}

        /** \brief Assignment operator for an InOutBlockData instance */
        InOutBlockData& operator=(const InOutBlockData& other)
        {
            this->m_idx = other.m_idx;
            return *this;
        }

    public:     // API for a BlockData

        /**
         * \brief Predicate to determine if the associated block is a leaf
         *
         * \return True, if the block is a leaf, False otherwise
         */
        bool isLeaf() const { return m_idx > INTERNAL_BLOCK; }

        /** \brief Marks the associated block as internal */
        void setInternal()  { m_idx = INTERNAL_BLOCK; }

        /** \brief Marks the associated block as a non-block (i.e. not in the tree) */
        void setNonBlock()  { m_idx = NON_BLOCK; }

        /**
         * \brief Predicate to determine if the associated block is in the tree
         *
         * \return True, if the block is in the tree (internal or leaf), False otherwise
         */
        bool isBlock() const { return m_idx != NON_BLOCK; }

    public:     // Other functions

        /**
         * Clears the data associated with the block
         * \note This function is currently a no-op
         * */
        void clear()
        {
            // No-op for now -- eventually, will need to do something about the index
        }

        /**
         * Predicate to determine if the associated block has data (i.e. it is a gray block)
         * \return True, if the block has data, False otherwise
         * */
        bool hasData() const { return m_idx >= 0; }

        /**
         * Returns the index of the data associated with the block
         */
        const int& dataIndex() const
        {
            //SLIC_ASSERT(hasData());
            return m_idx;
        }

        /**
         * \brief Sets the block as gray, and provides the index of its associated data
         *
         * \param idx The index of the data associated with the gray leaf block
         * \pre The block must be a leaf block
         * \pre The passed in index, idx, must be a non-negative integer
         */
        void setGray(int idx)
        {
            SLIC_ASSERT(isLeaf());
            SLIC_ASSERT(idx >= 0);
            m_idx = idx;
        }

        /** Marks the block as Black (the entire domain is inside the surface) */
        void setBlack()
        {
            SLIC_ASSERT( isLeaf() );
            m_idx = LEAF_BLOCK_BLACK;
        }

        /** Marks the block as Black (the entire domain is outside the surface) */
        void setWhite()
        {
            SLIC_ASSERT( isLeaf() );
            m_idx = LEAF_BLOCK_WHITE;
        }

        /** Sets the data associated with the block to the given index idx */
        void setData(int idx)
        {
            m_idx = idx;
        }

        /** Marks the block as uncolored */
        void setUncoloredLeaf()
        {
            SLIC_ASSERT( isLeaf() );
            m_idx = LEAF_BLOCK_UNCOLORED;
        }

        /**
         * \brief Find the 'color' of this LeafBlock
         *
         * 'Black' indicates that the entire block is within the surface
         * 'White' indicates that the entire block is outside the surface
         * 'Gray' indicates that the block intersects the surface geometry
         * Leaves that haven't been colored yet are 'Undetermined'
         */
        LeafColor color() const
        {
            if(hasData())
                return Gray;

            switch( m_idx )
            {
            case LEAF_BLOCK_BLACK:     return Black;
            case LEAF_BLOCK_WHITE:     return White;
            case LEAF_BLOCK_UNCOLORED: return Undetermined;
            }

            SLIC_ASSERT_MSG(false, "Invalid state in InOuLeafData::color()");
            return Undetermined;
        }

        /** Predicate to determine if the associated block has a color
         * \return True if the block has a color, false otherwise
         * \sa color()
         */
        bool isColored() const
        {
            return color() != Undetermined;
        }

        /** Friend function to compare equality of two InOutBlockData instances  */
        friend bool operator==(const InOutBlockData& lhs, const InOutBlockData& rhs)
        {
            return lhs.m_idx == rhs.m_idx;
        }
    private:
        int m_idx;
    };

    /**
     * Free function to print an InOutBlockData to an output stream
     * \param os The output stream to write to
     * \param iob The InOUtBlockData instance that we are writing
     */
    std::ostream& operator<<(std::ostream& os, const InOutBlockData& iob)
    {
        os << "InOutBlockData{"
           << "isLeaf: " << (iob.isLeaf() ? "yes" : "no");

        bool showData = true;

        if(iob.isLeaf())
        {
            os <<", color: ";
            switch(iob.color())
            {
            case InOutBlockData::Gray:
                os <<"Gray";
                break;
            case InOutBlockData::White:
                os <<"White";
                showData = false;
                break;
            case InOutBlockData::Black:
                os <<"Black";
                showData = false;
                break;
            default:
                os <<"Undetermined";
                break;
            }
        }

        if(showData)
        {
            os<<", dataIndex: ";
            if(!iob.hasData())
                os << "<no data>";
            else
                os << iob.dataIndex();
        }

        os << "}";

        return os;
    }

    /**
     * \brief Verbose BlockDataType for an InOutOctree
     *
     * \note Used when generating the octree.
     */
    class DynamicGrayBlockData
    {
    public:
        enum { NO_VERTEX = -1 };

        typedef int VertexIndex;
        typedef int TriangleIndex;

        typedef std::vector<TriangleIndex> TriangleList;

    public:
        /**
         * \brief Default constructor for an InOutLeafData
         */
        DynamicGrayBlockData() : m_vertIndex(NO_VERTEX), m_isLeaf(true) {}

        /**
         * \brief Constructor for an InOutLeafData
         *
         * \param vInd The index of a vertex (optional; default is to not set a vertex)
         */
        DynamicGrayBlockData(VertexIndex vInd, bool isLeaf): m_vertIndex(vInd), m_isLeaf(isLeaf) {}

        /**
         * \brief Copy constructor for an DynamicGrayBlockData instance
         */
        DynamicGrayBlockData(const DynamicGrayBlockData& other)
            : m_vertIndex(other.m_vertIndex)
            , m_tris(other.m_tris)
            , m_isLeaf(other.m_isLeaf)
        {}

        /**
         * \brief Assignment operator for an InOutLeafData instance
         */
        DynamicGrayBlockData& operator=(const DynamicGrayBlockData& other)
        {
            this->m_vertIndex = other.m_vertIndex;

            this->m_tris.reserve( other.m_tris.size());
            std::copy(other.m_tris.begin(), other.m_tris.end(), std::back_inserter(this->m_tris));

            this->m_isLeaf = other.m_isLeaf;

            return *this;
        }

//        /**
//         * \brief Removes all indexed data from this leaf
//         */
//        void clear()
//        {
//            m_isLeaf = false;
//            m_vertIndex = NO_VERTEX;
//            m_tris.clear();
//            m_tris = TriangleList(0);    // reconstruct to deallocate memory
//        }

        /**
         * \brief Equality operator to determine if two DynamicGrayBlockData instances are equivalent
         */
        friend bool operator==(const DynamicGrayBlockData& lhs, const DynamicGrayBlockData& rhs )
        {
            return //(static_cast<const BlockData&>(lhs) == static_cast<const BlockData&>(rhs))
                //&&
                (lhs.m_vertIndex == rhs.m_vertIndex)
                && (lhs.m_tris.size() == rhs.m_tris.size())     // Note: We are not checking the contents
                // && (lhs.m_tris == rhs.m_tris)                //       of the triangle array, only the size
                && lhs.m_isLeaf == rhs.m_isLeaf
                ;
        }


    public: // Functions related to whether this is a leaf

        /** Predicate to determine if the associated block is a leaf in the octree */
        bool isLeaf() const { return m_isLeaf; }

        /** Sets a flag to determine whether the associated block is a leaf or internal */
        void setLeafFlag(bool isLeaf) { m_isLeaf = isLeaf; }

    public: // Functions related to the associated vertex

        /**
         * \brief Checks whether there is a vertex associated with this leaf
         */
        bool hasVertex() const { return m_vertIndex >= 0; }

        /** Sets the vertex associated with this leaf */
        void setVertex(VertexIndex vInd) { m_vertIndex = vInd; }

        /** Clears the associated vertex index */
        void clearVertex() { m_vertIndex = NO_VERTEX; }

        /** Accessor for the index of the vertex associated with this leaf */
        VertexIndex& vertexIndex() { return m_vertIndex; }

        /** Constant accessor for the index of the vertex associated with this leaf */
        const VertexIndex& vertexIndex() const { return m_vertIndex; }


    public: // Functions related to the associated triangles
        /** Check whether this Leaf has any associated triangles */
        bool hasTriangles() const { return !m_tris.empty(); }

        /**
         * Reserves space for a given number of triangles
         * \param count The number of triangles for which to reserve space
         */
        void reserveTriangles(int count) { m_tris.reserve(count); }

        /** Find the number of triangles associated with this leaf */
        int numTriangles() const { return m_tris.size(); }

        /** Associates the surface triangle with the given index with this block */
        void addTriangle(TriangleIndex tInd) { m_tris.push_back(tInd); }

        /** Returns a const reference to the list of triangle indexes associated with the block */
        const TriangleList& triangles() const { return m_tris; }

        /** Returns a reference to the list of triangle indexes associated with the block */
        TriangleList& triangles() { return m_tris;}

    private:
        VertexIndex m_vertIndex;
        TriangleList m_tris;
        bool m_isLeaf;
    };

    /**
     * Free function to print a DynamicGrayBlockData instance to an output stream
     */
    std::ostream& operator<<(std::ostream& os, const DynamicGrayBlockData& bData)
    {
        os << "DynamicGrayBlockData{";

        os <<"isLeaf: " << ( bData.isLeaf() ?  "yes" : "no");

        os <<", vertex: ";
        if( bData.hasVertex() )
            os << bData.vertexIndex();
        else
            os << "<none>";

        os <<", triangles: ";
        if( bData.hasTriangles() )
        {
            int numTri = bData.numTriangles();
            os << "(" <<  numTri << ") {";
            for(int i=0; i< numTri; ++i)
                os << bData.triangles()[i] << ( (i == numTri-1) ? "} " : ",");
        }

        os<< "}";

        return os;
    }


/**
 * \class
 * \brief Handles generation of a point containment spatial index over a surface mesh
 *
 * The point containment queries determines whether a given arbitrary point in space
 * lies inside or outside of the surface
 */
template<int DIM>
class InOutOctree : public SpatialOctree<DIM, InOutBlockData>
{
public:

    typedef OctreeBase<DIM, InOutBlockData> OctreeBaseType;
    typedef SpatialOctree<DIM, InOutBlockData> SpatialOctreeType;

    typedef typename SpatialOctreeType::GeometricBoundingBox GeometricBoundingBox;
    typedef typename SpatialOctreeType::SpacePt SpacePt;
    typedef typename SpatialOctreeType::SpaceVector SpaceVector;
    typedef typename SpatialOctreeType::BlockIndex BlockIndex;
    typedef typename OctreeBaseType::GridPt GridPt;


private:
    enum GenerationState { INOUTOCTREE_UNINITIALIZED
                         , INOUTOCTREE_VERTICES_INSERTED
                         , INOUTOCTREE_MESH_REORDERED
                         , INOUTOCTREE_ELEMENTS_INSERTED
                         , INOUTOCTREE_LEAVES_COLORED
    };

    /**
     * \brief A utility class that wraps the access to the mesh data
     *
     * This class helps separate the specifics of accessing the underlying mesh for on InOutOctree.
     * It is customized for unstructured Triangle meshes, but we will later want to apply
     * the InOutOctree to other mesh types, e.g. Segments and Nurbs in 2D, Bilinear quads and NURBS in 3D.
     *
     * We can later specialize this class for triangle meshes and implement other customized mesh wrappers.
     */
    class MeshWrapper
    {
    public:
        typedef int VertexIndex;
        typedef int TriangleIndex;
        typedef axom::mint::Mesh SurfaceMesh;

        /** \brief A vertex index to indicate that there is no associated vertex */
        static const VertexIndex NO_VERTEX = -1;

        /** \brief A constant for the number of boundary vertices in a triangle */
        static const int NUM_TRI_VERTS = 3;

        typedef axom::primal::Triangle<double, DIM> SpaceTriangle;


        typedef axom::slam::PositionSet MeshVertexSet;
        typedef axom::slam::PositionSet MeshElementSet;

        typedef axom::slam::Map<VertexIndex> VertexIndexMap;
        typedef axom::slam::Map<SpacePt> VertexPositionMap;

        typedef axom::slam::policies::CompileTimeStrideHolder<VertexIndex, NUM_TRI_VERTS>  TVStride;
        typedef axom::slam::StaticConstantRelation<TVStride, MeshElementSet, MeshVertexSet> TriangleVertexRelation;
        typedef typename TriangleVertexRelation::RelationSet TriVertIndices;

    public:
        /** \brief Constructor for a mesh wrapper */
        MeshWrapper(SurfaceMesh*& meshPtr)
            : m_surfaceMesh(meshPtr)
            , m_vertexSet(0)
            , m_elementSet(0)
            , m_vertexPositions(&m_vertexSet)
            , m_triangleToVertexRelation()
            , m_meshWasReindexed(false)
        {}

        /** Const accessor to the vertex set of the wrapped surface mesh */
        const MeshVertexSet& vertexSet() const { return m_vertexSet; }

        /** Const accessor to the element set of the wrapped surface mesh */
        const MeshElementSet& elementSet() const { return m_elementSet; }

        /** Accessor to the vertex set of the wrapped surface mesh */
        MeshVertexSet& vertexSet() { return m_vertexSet; }

        /** Accessor to the element set of the wrapped surface mesh */
        MeshElementSet& elementSet() { return m_elementSet; }

        /** Accessor for the number of vertices in the wrapped surface mesh */
        int numMeshVertices() const
        {
            if(m_meshWasReindexed)
                return m_vertexSet.size();
            else
                return m_surfaceMesh->getMeshNumberOfNodes();
        }

        /** Accessor for the number of elements in the wrapped surface mesh */
        int numMeshElements() const
        {
            if(m_meshWasReindexed)
                return m_elementSet.size();
            else
                return m_surfaceMesh->getMeshNumberOfCells();
        }

        /** Predicate to determine if the wrapped surface mesh has been reindexed */
        bool meshWasReindexed() const { return m_meshWasReindexed;}


        /**
         * \brief Helper function to retrieve the position of the vertex from the mesh
         * \param idx The index of the vertex within the surface mesh
         */
        SpacePt getMeshVertexPosition(VertexIndex idx) const
        {
            if(m_meshWasReindexed)
            {
                return m_vertexPositions[idx];
            }
            else
            {
                SpacePt pt;
                m_surfaceMesh->getMeshNode(idx, pt.data() );
                return pt;
            }
        }

        /** Returns the spatial position of the vertex with index idx of the wrapped surface mesh */
        const SpacePt& vertexPosition(VertexIndex idx) const
        {
            return m_vertexPositions[idx];
        }


        /**
         * Returns the indices of the boundary vertices of the element of the wrapped surface mesh with the given index
         * \param idx The index of an element within the surface mesh
         */
        TriVertIndices triangleVertexIndices(TriangleIndex idx) const
        {
            return m_triangleToVertexRelation[idx];
        }

        /**
         * \brief Helper function to compute the bounding box of a triangle
         *
         * \param idx The triangle's index within the surface mesh
         */
        GeometricBoundingBox triangleBoundingBox(TriangleIndex idx) const
        {
            // Get the ids of the verts bounding this triangle
            TriVertIndices vertIds = triangleVertexIndices(idx);

            GeometricBoundingBox bb(vertexPosition( vertIds[0] ));
            bb.addPoint(vertexPosition( vertIds[1] ));
            bb.addPoint(vertexPosition( vertIds[2] ));

            return bb;
        }

        /**
         * \brief Utility function to retrieve the positions of the triangle's vertices
         *
         * \return A triangle instance whose vertices are positioned in space
         */
        SpaceTriangle trianglePositions(TriangleIndex idx) const
        {
            TriVertIndices verts = triangleVertexIndices(idx);
            return SpaceTriangle( vertexPosition(verts[0])
                                , vertexPosition(verts[1])
                                , vertexPosition(verts[2]));
        }

        /**
         * \brief Checks whether the indexed triangle contains a reference to the given vertex
         */
        bool incidentInVertex(const TriVertIndices& triVerts, VertexIndex vIdx) const
        {
            return (triVerts[0] == vIdx) || (triVerts[1] == vIdx) || (triVerts[2] == vIdx);
        }


        /**
         * \brief Finds the index of a vertex in triangle t1 that is not in triangle t0
         *
         * \param t0 The index of the first triangle
         * \param t1 The index of the second triangle
         * \pre t0 and t1 must be distinct triangles
         * \return The index of a vertex in t1 that is not in t0,
         *         (NO_VERTEX if one does not exist)
         */
        VertexIndex distinctVertex(TriangleIndex t0, TriangleIndex t1) const
        {
            SLIC_ASSERT_MSG(t0 != t1, "Expected two different triangle indices in "
                            <<"quest::MeshWrapper::findDistinctVertex, got triangle index " << t0 << " twice.");

            // Find a vertex from the local surface that is not incident in the first triangle
            TriVertIndices tvRel0 = triangleVertexIndices( t0 );
            TriVertIndices tvRel1 = triangleVertexIndices( t1 );

            for(int i=0; i< NUM_TRI_VERTS; ++i)
                if(!incidentInVertex(tvRel0, tvRel1[i]))
                    return tvRel1[i];

            SLIC_ASSERT_MSG(false, "There should be a vertex in triangle "
                                << t1 << " that was not in " << t0);
            return NO_VERTEX;
        }

        /**
         * \brief Determine if the two given triangles have a vertex in common
         *
         * \param t0 The index of the first triangle
         * \param t1 The index of the second triangle
         * \param [out] sharedVert The index of the shared vertex, if it exists
         * \return true if the two triangles have a vertex in common (returned in sharedVert), false otherwise
         */
        bool haveSharedVertex(TriangleIndex t0, TriangleIndex t1, VertexIndex& sharedVert) const
        {
            // There are two triangles -- check that they have at least one common vertex
            TriVertIndices tvRel0 = triangleVertexIndices(t0);
            TriVertIndices tvRel1 = triangleVertexIndices(t1);

            for(int i=0; i< NUM_TRI_VERTS; ++i)
            {
                if(incidentInVertex(tvRel0, tvRel1[i]))
                {
                    sharedVert = tvRel1[i];
                    return true;
                }
            }
            return false;
        }

        /**
         * \brief Determine if the three given triangles have a vertex in common
         *
         * \param t0 The index of the first triangle
         * \param t1 The index of the second triangle
         * \param t2 The index of the second triangle
         * \param [out] sharedVert The index of the shared vertex, if it exists
         * \return true if the three triangles have a vertex in common (returned in sharedVert), false otherwise
         */
        bool haveSharedVertex(TriangleIndex t0, TriangleIndex t1, TriangleIndex t2, VertexIndex& sharedVert) const
        {
            TriVertIndices t0Verts = triangleVertexIndices(t0);
            TriVertIndices t1Verts = triangleVertexIndices(t1);
            TriVertIndices t2Verts = triangleVertexIndices(t2);

            for(int i=0; i< NUM_TRI_VERTS; ++i)
            {
                // check if a vertex from the third triangle is in the first and second
                if(incidentInVertex( t0Verts, t2Verts[i]) && incidentInVertex( t1Verts, t2Verts[i]) )
                {
                    sharedVert = t2Verts[i];
                    return true;
                }
            }

            return false;
        }


        /**
         * \brief Reindexes the mesh vertices and triangle indices using the given map
         *
         * \param numVertices The number of vertices in the new mesh
         * \param vertexIndexMap A mapping from the old vertex indices to the new ones
         * \note This step clears out the original mesh,
         * which can be reconstructed using the regenerateSurfaceMesh() function
         */
        void reindexMesh(int numVertices, const VertexIndexMap& vertexIndexMap)
        {
            // Create a vertex set on the new vertices and grab coordinates from the old ones
            m_vertexSet = MeshVertexSet( numVertices );
            m_vertexPositions = VertexPositionMap(&m_vertexSet);

            int numOrigVertices = numMeshVertices();
            for(int i=0; i< numOrigVertices; ++i)
            {
                const VertexIndex& vInd = vertexIndexMap[i];
                m_vertexPositions[vInd] = getMeshVertexPosition(i);
            }

            // Update the vertex IDs of the triangles to the new vertices and create a SLAM relation on these
            int numOrigTris = numMeshElements();
            std::vector<int> tvrefs;
            tvrefs.reserve(NUM_TRI_VERTS * numOrigTris);
            for(int i=0; i< numOrigTris ; ++i)
            {
                // Grab relation from mesh
                int vertIds[NUM_TRI_VERTS];
                m_surfaceMesh->getMeshCell(i, vertIds);

               // Remap the vertex IDs
               for(int j=0; j< NUM_TRI_VERTS; ++j)
                   vertIds[j] = vertexIndexMap[ vertIds[j] ];

               // Add to relation if not degenerate triangles (namely, we need 3 unique vertex IDs)
               if(    (vertIds[0] != vertIds[1])
                   && (vertIds[1] != vertIds[2])
                   && (vertIds[2] != vertIds[0]) )
               {
                   tvrefs.push_back (vertIds[0]);
                   tvrefs.push_back (vertIds[1]);
                   tvrefs.push_back (vertIds[2]);
               }
            }

            m_elementSet = MeshElementSet( tvrefs.size() / NUM_TRI_VERTS );
            m_triangleToVertexRelation = TriangleVertexRelation(&m_elementSet, &m_vertexSet);
            m_triangleToVertexRelation.bindRelationData(tvrefs);


            // Delete old mesh, and NULL its pointer
            delete m_surfaceMesh;
            m_surfaceMesh = AXOM_NULLPTR;

            m_meshWasReindexed = true;
        }

        /** \brief Add the vertex positions and triangle boundary relations to the surface mesh */
        void regenerateSurfaceMesh()
        {
          if(m_surfaceMesh != AXOM_NULLPTR)
          {
              delete m_surfaceMesh;
              m_surfaceMesh = AXOM_NULLPTR;
          }

          typedef axom::mint::UnstructuredMesh< MINT_TRIANGLE > TriangleMesh;
          TriangleMesh* triMesh = new TriangleMesh(3);

          // Add vertices to the mesh (i.e. vertex positions)
          for(int i=0; i< m_vertexSet.size(); ++i)
          {
              const SpacePt& pt = vertexPosition(i);
              triMesh->insertNode(pt[0], pt[1], pt[2]);
          }

          // Add triangles to the mesh (i.e. boundary vertices)
          for(int i=0; i< m_elementSet.size(); ++i)
          {
              const TriangleIndex* tv = &triangleVertexIndices(i)[0];
              triMesh->insertCell(tv, MINT_TRIANGLE, NUM_TRI_VERTS);
          }

          m_surfaceMesh = triMesh;
        }

    private:
        SurfaceMesh*& m_surfaceMesh; /** pointer to pointer to allow changing the mesh */

        MeshVertexSet m_vertexSet;
        MeshElementSet m_elementSet;
        VertexPositionMap m_vertexPositions;
        TriangleVertexRelation m_triangleToVertexRelation;

        bool m_meshWasReindexed;
    };


public:

    typedef typename MeshWrapper::SurfaceMesh SurfaceMesh;
    typedef axom::mint::UnstructuredMesh< MINT_MIXED_CELL > DebugMesh;

    typedef typename MeshWrapper::VertexIndex VertexIndex;
    typedef typename MeshWrapper::TriangleIndex TriangleIndex;
    typedef typename MeshWrapper::SpaceTriangle SpaceTriangle;

    typedef typename MeshWrapper::MeshVertexSet MeshVertexSet;
    typedef typename MeshWrapper::TriVertIndices TriVertIndices;
    typedef typename MeshWrapper::VertexIndexMap VertexIndexMap;



    // Some type defs for the Relations from Gray leaf blocks to mesh vertices and elements
    static const int MAX_VERTS_PER_BLOCK = 1;
    typedef axom::slam::Map<BlockIndex> VertexBlockMap;

    typedef axom::slam::PositionSet GrayLeafSet;
    typedef axom::slam::policies::CompileTimeStrideHolder<VertexIndex, MAX_VERTS_PER_BLOCK>  BVStride;
    typedef axom::slam::StaticConstantRelation<BVStride, GrayLeafSet, MeshVertexSet> GrayLeafVertexRelation;
    typedef axom::slam::StaticVariableRelation                                       GrayLeafElementRelation;
    typedef typename GrayLeafElementRelation::RelationSet TriangleIndexSet;

    typedef axom::slam::Map<GrayLeafSet>            GrayLeafsLevelMap;
    typedef axom::slam::Map<GrayLeafVertexRelation> GrayLeafVertexRelationLevelMap;
    typedef axom::slam::Map<GrayLeafElementRelation> GrayLeafElementRelationLevelMap;

public:
    /**
     * \brief Construct an InOutOctree to handle containment queries on a surface mesh
     *
     * \param [in] bb The spatial extent covered by the octree
     * \note We slightly scale the bounding box so all mesh elements are guaranteed
     *       to be enclosed by the octree
     */
    InOutOctree(const GeometricBoundingBox& bb, SurfaceMesh*& meshPtr)
        : SpatialOctreeType( GeometricBoundingBox(bb).scale(1.0001) )
        , m_meshWrapper(meshPtr)
        , m_vertexToBlockMap(&m_meshWrapper.vertexSet())
        //
        , m_grayLeafsMap( &this->m_levels)
        , m_grayLeafToVertexRelationLevelMap( &this->m_levels )
        , m_grayLeafToElementRelationLevelMap( &this->m_levels )
        //
        , m_generationState( INOUTOCTREE_UNINITIALIZED)
    {
    }

    /**
     * \brief Generate the spatial index over the triangle mesh
     */
    void generateIndex();



    /**
     * \brief The point containment query.
     *
     * \param pt The point that we want to check for containment within the surface
     * \return True if the point is within (or on) the surface, false otherwise
     * \note Points outside the octree's bounding box are considered outside the surface
     */
    bool within(const SpacePt& pt) const;

private:

    /**
     * \brief Helper function to insert a vertex into the octree
     *
     * \param idx The index of the vertex that we are inserting
     * \param startingLevel (optional, default = 0) The octree level at which
     * to begin the search for the leaf node covering this vertex
     */
    void insertVertex(VertexIndex idx, int startingLevel = 0);

    /**
     * \brief Insert all triangles of the mesh into the octree, generating a PM octree
     */
    void insertMeshTriangles();


    /**
     * \brief Set a color for each leaf block of the octree.
     *
     * Black blocks are entirely within the surface, white blocks are entirely outside the surface
     * and Gray blocks intersect the surface.
     */
    void colorOctreeLeaves();


    /**
     * Use octree index over mesh vertices to convert the 'triangle soup'
     * from the stl file into an indexed triangle mesh representation.
     * In particular, all vertices in the mesh that are nearly coincident will be merged,
     * and degenerate triangles (where the three vertices do not have unique indices)
     * will be removed.
     */
    void updateSurfaceMeshVertices();

private:

    /**
     * \brief Checks if all indexed triangles in the block share a common vertex
     *
     * \param leafBlock [in] The current octree block
     * \param leafData [inout] The data associated with this block
     * \note A side effect of this function is that we set the leafData's vertex to the common
     * vertex if one is found
     * \return True, if all triangles indexed by this leaf share a common vertex, false otherwise.
     */
    bool allTrianglesIncidentInCommonVertex(const BlockIndex& leafBlock, DynamicGrayBlockData& leafData) const;

    /**
     * \brief Finds a color for the given block blk and propagates to neighbors
     *
     * \note Propagates color to same-level and coarser level neighbors
     * \param blk The block to color
     * \param blkData The data associated with this block
     * \return True if we were able to find a color for blk, false otherwise
     */
    bool colorLeafAndNeighbors(const BlockIndex& blk, InOutBlockData& blkData);


    /**
     * \brief Predicate to determine if the vertex is indexed by the blk
     *
     * \pre This function assumes the vertices have been inserted and the mesh has been reordered
     * \param vIdx The index of the vertex to check
     * \param blk The block that we are checking against
     * \return true if vIdx is indexed by blk, false otherwise
     */
    bool blockIndexesVertex(VertexIndex vIdx, const BlockIndex& blk) const
    {
        SLIC_ASSERT(m_generationState >= INOUTOCTREE_MESH_REORDERED);
        return vIdx >=0 && m_vertexToBlockMap[vIdx] == blk;
    }

    /**
     * \brief Predicate to determine if any of the elements vertices are indexed by the given BlockIndex
     *
     * \pre This function assumes the vertices have been inserted and the mesh has been reordered
     * \param tIdx The index of the triangle to check
     * \param blk The block that we are checking against
     * \return true if one of the triangle's vertices are indexed by blk, false otherwise
     */
    bool blockIndexesElementVertex(TriangleIndex tIdx, const BlockIndex& blk) const
    {
        SLIC_ASSERT(m_generationState >= INOUTOCTREE_MESH_REORDERED);

        TriVertIndices tVerts = m_meshWrapper.triangleVertexIndices(tIdx);
        for(int i=0; i< tVerts.size(); ++i)
        {
            // Note: Use the vertex to block cache to avoid numerical degeneracies
            if( blk == m_vertexToBlockMap[ tVerts[i]])
                return true;
        }
        return false;
    }

    /**
     * \brief Determines whether the specified point is within the gray leaf
     *
     * \param pt The point we are querying
     * \param leafBlk The block of the gray leaf
     * \param data The data associated with the leaf block
     * \return True, if the point is inside the local surface associated with this block, false otherwise
     */
    bool withinGrayBlock( const SpacePt & pt, const BlockIndex& leafBlk, const InOutBlockData&  data) const;

    /**
     * \brief Returns the index of the mesh vertex associated with the given leaf block
     *
     * \pre leafBlk is a leaf block of the octree
     * \param leafBlk The BlockIndex of a leaf block in the octree
     * \param leafData The data associated with this leaf block
     * \return The index of the mesh vertex associated with this leaf block
     */
    VertexIndex leafVertex(const BlockIndex& leafBlk, const InOutBlockData&  leafData) const;

    /**
     * \brief Returns the set of mesh triangle indices associated with the given leaf block
     *
     * \pre leafBlk is a leaf block of the octree
     * \param leafBlk The BlockIndex of a leaf block in the octree
     * \param leafData The data associated with this leaf block
     * \return The set of mesh triangle indices associated with this leaf block
     */
    TriangleIndexSet leafTriangles(const BlockIndex& leafBlk, const InOutBlockData&  leafData) const;


private:
    DISABLE_COPY_AND_ASSIGNMENT(InOutOctree);
    DISABLE_MOVE_AND_ASSIGNMENT(InOutOctree);

    /** \brief Checks internal consistency of the octree representation */
    void checkValid() const;


    /** \brief Helper function to verify that all leaves at the given level have a color */
    void checkAllLeavesColoredAtLevel(int AXOM_DEBUG_PARAM(level)) const;

    void dumpOctreeMeshVTK( const std::string& name) const;


    /**
     * \brief Utility function to dump a single element vtk file
     */
    void dumpMeshVTK( const std::string& name
                    , int idx
                    , const BlockIndex& block
                    , const GeometricBoundingBox& blockBB
                    , bool isTri
                    ) const;

    /**
     * \brief Utility function to dump any Inside blocks whose neighbors are outside (and vice-versa)
     *
     * \note There should not be any such blocks in a valid InOutOctree
     */
    void dumpDifferentColoredNeighborsMeshVTK(const std::string& name) const;



    /**
     * Utility function to print some statistics about the InOutOctree instance
     */
    void printOctreeStats() const;


protected:

    MeshWrapper m_meshWrapper;

    VertexBlockMap m_vertexToBlockMap;

    GrayLeafsLevelMap m_grayLeafsMap;
    GrayLeafVertexRelationLevelMap m_grayLeafToVertexRelationLevelMap;
    GrayLeafElementRelationLevelMap m_grayLeafToElementRelationLevelMap;

    GenerationState m_generationState;
};


namespace{
#ifdef AXOM_DEBUG
    /**
     * \brief Utility function to print the vertex indices of a cell
     */
    std::ostream& operator<<(std::ostream& os, const InOutOctree<3>::TriVertIndices& tvInd)
    {
        os<<"[";
        for(int i=0; i< tvInd.size(); ++i)
            os << tvInd[i] << ( (i == tvInd.size()-1) ? "]" : ",");

        return os;
    }

    std::ostream& operator<<(std::ostream& os, const InOutOctree<3>::TriangleIndexSet& tSet)
    {
        os<<"[";
        for(int i=0; i< tSet.size(); ++i)
            os << tSet[i] << ( (i == tSet.size()-1) ? "]" : ",");

        return os;
    }

#endif
}

template<int DIM>
void InOutOctree<DIM>::generateIndex ()
{
    typedef axom::utilities::Timer Timer;

    // Loop through mesh vertices
    SLIC_INFO("  Generating InOutOctree over surface mesh with "
              << m_meshWrapper.numMeshVertices()  << " vertices and "
              << m_meshWrapper.numMeshElements() << " elements.");

    Timer timer;

    // STEP 1 -- Add mesh vertices to octree
    timer.start();
    int numMeshVerts = m_meshWrapper.numMeshVertices();
    for(int idx=0; idx < numMeshVerts; ++idx)
    {
        insertVertex(idx);
    }
    timer.stop();
    m_generationState = INOUTOCTREE_VERTICES_INSERTED;
    SLIC_INFO("\t--Inserting vertices took " << timer.elapsed() << " seconds.");



    // STEP 1(b) -- Update the mesh vertices and cells with after vertex welding from octree
    timer.start();
    updateSurfaceMeshVertices();
    timer.stop();
    m_generationState = INOUTOCTREE_MESH_REORDERED;

    SLIC_INFO("\t--Updating mesh took " << timer.elapsed() << " seconds.");
    SLIC_INFO("  After inserting vertices, reindexed mesh has "
              << m_meshWrapper.numMeshVertices() << " vertices" << " and "
              << m_meshWrapper.numMeshElements() << " triangles." );


  #ifdef DUMP_OCTREE_INFO
    // -- Print some stats about the octree
    SLIC_INFO("** Octree stats after inserting vertices");
    dumpOctreeMeshVTK("prOctree.vtk");
    printOctreeStats();
  #endif
    checkValid();


    // STEP 2 -- Add mesh triangles to octree
    timer.start();
    insertMeshTriangles();
    timer.stop();
    m_generationState = INOUTOCTREE_ELEMENTS_INSERTED;
    SLIC_INFO("\t--Inserting triangles took " << timer.elapsed() << " seconds.");


    // STEP 3 -- Color the blocks of the octree -- Black (in), White(out), Gray(Intersects surface)
    timer.start();
    colorOctreeLeaves();

    timer.stop();
    m_generationState = INOUTOCTREE_LEAVES_COLORED;
    SLIC_INFO("\t--Coloring octree leaves took " << timer.elapsed() << " seconds.");

    // -- Print some stats about the octree
  #ifdef DUMP_OCTREE_INFO
    SLIC_INFO("** Octree stats after inserting triangles");
    dumpOctreeMeshVTK("pmOctree.vtk");
    dumpDifferentColoredNeighborsMeshVTK("differentNeighbors.vtk");
    printOctreeStats();
  #endif
    checkValid();


    // CLEANUP -- Finally, fix up the surface mesh after octree operations
    timer.start();
    m_meshWrapper.regenerateSurfaceMesh();
    timer.stop();
    SLIC_INFO("\t--Regenerating the mesh took " << timer.elapsed() << " seconds.");

    SLIC_INFO("  Finished generating the InOutOctree.");
}


template<int DIM>
void InOutOctree<DIM>::insertVertex (VertexIndex idx, int startingLevel)
{
    static const double EPS_SQ = 1e-18;


    const SpacePt pt = m_meshWrapper.getMeshVertexPosition(idx);

    BlockIndex block = this->findLeafBlock(pt, startingLevel);
    InOutBlockData& blkData = (*this)[block];


    QUEST_OCTREE_DEBUG_LOG_IF( idx == DEBUG_VERT_IDX
              , "\t -- inserting pt " << pt << " with index " << idx
                  << ". Looking at block " << block
                  << " w/ blockBB " << this->blockBoundingBox(block)
                  << " indexing leaf vertex " << blkData.dataIndex()
              );

    if(! blkData.hasData() )
    {
        blkData.setData(idx);

        // Update the vertex-to-block map for this vertex
        if(m_generationState >= INOUTOCTREE_MESH_REORDERED)
            m_vertexToBlockMap[idx] = block;
    }
    else
    {
        // check if we should merge the vertices
        VertexIndex origVertInd = blkData.dataIndex();
        if( squared_distance( pt, m_meshWrapper.getMeshVertexPosition(origVertInd) ) >= EPS_SQ )
        {
            blkData.clear();
            this->refineLeaf(block);

            insertVertex(origVertInd, block.childLevel() );
            insertVertex(idx, block.childLevel() );
        }
    }

    QUEST_OCTREE_DEBUG_LOG_IF( blkData.dataIndex() == DEBUG_VERT_IDX
              , "-- vertex " << idx << " is indexed in block " << block
                << ". Leaf vertex is " << blkData.dataIndex() );

}


template<int DIM>
void InOutOctree<DIM>::insertMeshTriangles ()
{
    typedef axom::utilities::Timer Timer;

    typedef typename OctreeBaseType::OctreeLevelType LeavesLevelMap;
    typedef typename OctreeBaseType::LevelMapIterator LeavesIterator;

    SLIC_ASSERT( m_meshWrapper.meshWasReindexed() );

    // Temporary arrays of DyamicGrayBlockData for current and next level
    typedef std::vector<DynamicGrayBlockData> DynamicLevelData;
    const int NUM_INIT_DATA_ENTRIES = 1<<10;
    DynamicLevelData currentLevelData;
    DynamicLevelData nextLevelData;
    currentLevelData.reserve( NUM_INIT_DATA_ENTRIES );
    nextLevelData.reserve( NUM_INIT_DATA_ENTRIES );

    // Raw data for relations from gray octree leaves to mesh vertices and elements
    std::vector<VertexIndex>   gvRelData;
    std::vector<int>           geSizeRelData;
    std::vector<TriangleIndex> geIndRelData;


    /// --- Initialize root level data
    BlockIndex rootBlock = this->root();
    InOutBlockData& rootData = (*this)[rootBlock];

    currentLevelData.push_back(DynamicGrayBlockData());
    DynamicGrayBlockData& dynamicRootData = currentLevelData[0];
    if( rootData.hasData() )
        dynamicRootData.setVertex( rootData.dataIndex() );
    dynamicRootData.setLeafFlag(rootData.isLeaf());
    rootData.setData(0);

    // Add all triangle references to the root
    int const numTris = m_meshWrapper.numMeshElements();
    dynamicRootData.triangles().reserve( numTris );
    for(int idx=0; idx < numTris; ++idx)
        dynamicRootData.addTriangle(idx);


    /// --- Iterate through octree levels and insert triangles into the blocks that they intersect
    for(int lev=0; lev< this->m_levels.size(); ++lev)
    {
        Timer levelTimer(true);

        gvRelData.clear();
        geSizeRelData.clear();
        geSizeRelData.push_back(0);
        geIndRelData.clear();

        int nextLevelDataBlockCounter = 0;

        LeavesLevelMap& levelLeafMap = this->getOctreeLevel( lev );
        for(LeavesIterator it=levelLeafMap.begin(), itEnd = levelLeafMap.end(); it != itEnd; ++it)
        {
            InOutBlockData& blkData = *it;

            if(! blkData.hasData() )
                continue;

            BlockIndex blk(it.pt(), lev);
            DynamicGrayBlockData& dynamicLeafData = currentLevelData[ blkData.dataIndex() ];

            bool isInternal = !dynamicLeafData.isLeaf();
            bool isLeafThatMustRefine = !isInternal && ! allTrianglesIncidentInCommonVertex(blk, dynamicLeafData);

            QUEST_OCTREE_DEBUG_LOG_IF(DEBUG_BLOCK_1 == blk || DEBUG_BLOCK_2 == blk
                 , "Attempting to insert triangles from block" << blk << "."
                 << "\n\tDynamic data: " << dynamicLeafData
                 << "\n\tBlock data: " << blkData
                 << "\n\tAbout to finalize? " << (!isInternal && !isLeafThatMustRefine  ? " yes" : "no" )
                );


            // Leaf blocks that don't refine are 'finalized' -- add  them to the current level's relations
            if( !isInternal && !isLeafThatMustRefine )
            {
                if( dynamicLeafData.hasTriangles())
                {
                    // Set the leaf data in the octree
                    blkData.setData( gvRelData.size() );

                    // Add the vertex index to the gray blocks vertex relation
                    gvRelData.push_back( dynamicLeafData.vertexIndex() );

                    // Add the triangles to the gray block's element relations
                    std::copy( dynamicLeafData.triangles().begin()
                             , dynamicLeafData.triangles().end()
                             , std::back_inserter(geIndRelData));
                    geSizeRelData.push_back( geIndRelData.size());

                    QUEST_OCTREE_DEBUG_LOG_IF(DEBUG_BLOCK_1 == blk || DEBUG_BLOCK_2 == blk
                         , "[Added block" << blk << " into tree as a gray leaf]."
                         << "\n\tDynamic data: " << dynamicLeafData
                         << "\n\tBlock data: " << blkData );
                }
            }
            else
            {
                /// Otherwise, we must distribute the block data among the children

                // Refine the leaf if necessary
                if( isLeafThatMustRefine )
                {
                    const VertexIndex vIdx = dynamicLeafData.vertexIndex();

                    this->refineLeaf(blk);
                    dynamicLeafData.setLeafFlag(false);

                    // Reinsert the vertex into the tree, if vIdx was indexed by blk
                    if( blockIndexesVertex(vIdx, blk) )
                        insertVertex(vIdx, blk.childLevel());
                }
                else if(isInternal)
                {
                    // Need to mark the leaf as internal since we were using its data
                    // as an index into the DynamicGrayBlockData array
                    blkData.setInternal();
                }

                SLIC_ASSERT_MSG( this->isInternal(blk)
                               , "Block " << blk <<" was refined, so it should be marked as internal." );


                /// Setup caches for data associated with children
                BlockIndex              childBlk[BlockIndex::NUM_CHILDREN];
                GeometricBoundingBox    childBB[BlockIndex::NUM_CHILDREN];
                DynamicGrayBlockData    childData[BlockIndex::NUM_CHILDREN];
                DynamicGrayBlockData*   childDataPtr[BlockIndex::NUM_CHILDREN];


                const typename LeavesLevelMap::BroodData& broodData = this->getOctreeLevel( lev+1 ).getBroodData(blk.pt());

                for(int j=0; j< BlockIndex::NUM_CHILDREN; ++j)
                {
                    childBlk[j]  = blk.child(j);
                    childBB[j]   = this->blockBoundingBox( childBlk[j] );

                    const InOutBlockData& childBlockData = broodData[j];
                    if(!childBlockData.hasData())
                    {
                        childData[j] = DynamicGrayBlockData();
                        childData[j].setLeafFlag( childBlockData.isLeaf() );
                    }
                    else
                    {
                        childData[j] = DynamicGrayBlockData(childBlockData.dataIndex(), childBlockData.isLeaf());
                    }

                    childDataPtr[j] = & childData[j];
                }

                // Check that the vector has enough capacity for all eight children
                // This ensures that our child data pointers will not be invalidated
                if(nextLevelData.capacity() < ( nextLevelData.size() + BlockIndex::NUM_CHILDREN ))
                    nextLevelData.reserve(nextLevelData.size() * 4);

                // Add all triangles to intersecting children blocks
                DynamicGrayBlockData::TriangleList& parentTriangles = dynamicLeafData.triangles();
                int numTriangles = parentTriangles.size();
                for(int i=0; i< numTriangles; ++i)
                {
                    TriangleIndex tIdx = parentTriangles[i];
                    SpaceTriangle spaceTri = m_meshWrapper.trianglePositions(tIdx);
                    GeometricBoundingBox tBB = m_meshWrapper.triangleBoundingBox(tIdx);

                    for(int j=0; j< BlockIndex::numChildren(); ++j)
                    {
                        bool shouldAddTriangle = childDataPtr[j]->isLeaf()
                                ? blockIndexesElementVertex(tIdx, childBlk[j]) || intersect(spaceTri, childBB[j])
                                : intersect(tBB, childBB[j])
                                ;

                        QUEST_OCTREE_DEBUG_LOG_IF( DEBUG_BLOCK_1 == childBlk[j]
                                                || DEBUG_BLOCK_2 == childBlk[j]  //&& tIdx == DEBUG_TRI_IDX
                             , "Attempting to insert triangle " << tIdx <<" @ " << spaceTri << " w/ BB " << tBB
                             << "\n\t into block" << childBlk[j] << " w/ BB " << childBB[j] << " and data " << *childDataPtr[j]
                             << "\n\tShould add? " << (shouldAddTriangle  ? " yes" : "no" )
                            );

                        if(shouldAddTriangle)
                        {
                            // Place the DynamicGrayBlockData in the array before adding its data
                            if(!childDataPtr[j]->hasTriangles())
                            {
                                // Copy the DynamicGrayBlockData into the array
                                nextLevelData.push_back(childData[j]);

                                // Update the child data pointer
                                childDataPtr[j] = &nextLevelData[nextLevelDataBlockCounter];

                                // Set the data in the octree to this index and update the index
                                (*this)[childBlk[j]].setData(nextLevelDataBlockCounter++);
                            }

                            childDataPtr[j]->addTriangle(tIdx);

                            QUEST_OCTREE_DEBUG_LOG_IF(DEBUG_BLOCK_1 == childBlk[j]
                                                   || DEBUG_BLOCK_2 == childBlk[j] //&& tIdx == DEBUG_TRI_IDX
                                 ,"Added triangle " << tIdx
                                 <<" @ " << spaceTri
                                 << " with verts " << m_meshWrapper.triangleVertexIndices(tIdx)
                                 << "\n\tinto block " << childBlk[j]
                                 << "\n\twith data " << *(childDataPtr[j])
                                 );
                        }
                    }
                }
            }
        }

        if(! levelLeafMap.empty() )
        {
            // Create the relations from gray leaves to mesh vertices and elements
            m_grayLeafsMap[lev] = GrayLeafSet( gvRelData.size() );

            m_grayLeafToVertexRelationLevelMap[lev] = GrayLeafVertexRelation(&m_grayLeafsMap[lev], &m_meshWrapper.vertexSet());
            m_grayLeafToVertexRelationLevelMap[lev].bindRelationData(gvRelData);

            m_grayLeafToElementRelationLevelMap[lev] = GrayLeafElementRelation(&m_grayLeafsMap[lev], &m_meshWrapper.elementSet());
            m_grayLeafToElementRelationLevelMap[lev].bindRelationData(geSizeRelData, geIndRelData);
        }

        currentLevelData.clear();
        nextLevelData.swap(currentLevelData);

        if(! levelLeafMap.empty() )
          SLIC_DEBUG("\tInserting triangles into level " << lev <<" took " << levelTimer.elapsed() <<" seconds.");
    }
}


template<int DIM>
void InOutOctree<DIM>::colorOctreeLeaves()
{
    // Note (KW): Existence of leaf implies that either
    // * it is gray
    // * one of its siblings is gray
    // * one of its siblings has a gray descendant

    typedef axom::utilities::Timer Timer;
    typedef typename OctreeBaseType::OctreeLevelType LeavesLevelMap;
    typedef typename OctreeBaseType::LevelMapIterator LeavesIterator;
    typedef std::vector<GridPt> GridPtVec;
    GridPtVec uncoloredBlocks;

    // Bottom-up traversal of octree
    for(int lev=this->maxLeafLevel()-1; lev >= 0; --lev)
    {
        uncoloredBlocks.clear();
        Timer levelTimer(true);

        LeavesLevelMap& levelLeafMap = this->getOctreeLevel(lev);
        for(LeavesIterator it=levelLeafMap.begin(), itEnd = levelLeafMap.end(); it != itEnd; ++it)
        {
            if( ! it->isLeaf() )
                continue;

            BlockIndex leafBlk(it.pt(), lev);
            InOutBlockData& blockData = *it;
            if(! colorLeafAndNeighbors( leafBlk, blockData) )
                uncoloredBlocks.push_back( leafBlk.pt());
        }

        // Iterate through the uncolored blocks until all have a color
        // This should terminate since we know that one of its siblings (or their descendants) is gray
        while(! uncoloredBlocks.empty())
        {
            int prevCount = uncoloredBlocks.size();
            AXOM_DEBUG_VAR(prevCount);

            GridPtVec prevVec;
            prevVec.swap(uncoloredBlocks);
            for(typename GridPtVec::iterator it = prevVec.begin(), itEnd = prevVec.end(); it < itEnd; ++it)
            {
                BlockIndex leafBlk(*it, lev);
                if(! colorLeafAndNeighbors( leafBlk, (*this)[leafBlk] ) )
                    uncoloredBlocks.push_back( *it);
            }

            SLIC_ASSERT_MSG( static_cast<int>(uncoloredBlocks.size()) < prevCount
                , "Problem coloring leaf blocks at level " << lev
                << ". There are  " << uncoloredBlocks.size() << " blocks that are still not colored."
                << " First problem block is: " << BlockIndex( uncoloredBlocks[0], lev)
                );
        }

        if(! levelLeafMap.empty() )
        {
            checkAllLeavesColoredAtLevel(lev);
            SLIC_DEBUG("\tColoring level "<< lev << " took " << levelTimer.elapsed() << " seconds.");
        }
    }
}


template<int DIM>
bool InOutOctree<DIM>::colorLeafAndNeighbors(const BlockIndex& leafBlk, InOutBlockData& leafData)
{
    bool isColored = leafData.isColored();

    QUEST_OCTREE_DEBUG_LOG_IF(leafBlk == DEBUG_BLOCK_1 || leafBlk == DEBUG_BLOCK_2
                              , "Trying to color " << leafBlk << " with data: " << leafData );


    if( ! isColored )
    {
        // Leaf does not yet have a color... try to find its color from same-level face neighbors
        for(int i=0; !isColored && i< leafBlk.numFaceNeighbors(); ++i)
        {
            BlockIndex neighborBlk = leafBlk.faceNeighbor(i);
            if( this->isLeaf( neighborBlk ) )
            {
                const InOutBlockData& neighborData = (*this)[neighborBlk];

                QUEST_OCTREE_DEBUG_LOG_IF(
                        DEBUG_BLOCK_1 == neighborBlk || DEBUG_BLOCK_2 == neighborBlk
                         || DEBUG_BLOCK_1 == leafBlk || DEBUG_BLOCK_2 == leafBlk
                     ,"Spreading color to block " << leafBlk << " with data " <<  leafData
                     << " and bounding box " << this->blockBoundingBox(leafBlk)
                     << " -- midpoint " << this->blockBoundingBox(leafBlk).centroid()
                     //
                     << "\n\t\t from " << neighborBlk << " with data " << neighborData
                     << " and bounding box " << this->blockBoundingBox(neighborBlk)
                     << " -- midpoint " << this->blockBoundingBox(neighborBlk).centroid()
                    );

                switch(neighborData.color())
                {
                case InOutBlockData::Black:
                    leafData.setBlack();
                    break;
                case InOutBlockData::White:
                    leafData.setWhite();
                    break;
                case InOutBlockData::Gray:
                    if( withinGrayBlock(this->blockBoundingBox(leafBlk).centroid(),  neighborBlk, neighborData) )
                        leafData.setBlack();
                    else
                        leafData.setWhite();
                    break;
                case InOutBlockData::Undetermined:
                    break;
                }

                isColored = leafData.isColored();

                QUEST_OCTREE_DEBUG_LOG_IF(isColored && (DEBUG_BLOCK_1 == neighborBlk || DEBUG_BLOCK_2 == neighborBlk)
                     , "Leaf block was colored -- " << leafBlk << " now has data " << leafData);
            }
        }
    }

    // If the block has a color, try to color its face neighbors at the same or coarser resolution
    if( isColored )
    {
        for(int i=0; i< leafBlk.numFaceNeighbors(); ++i)
        {
            BlockIndex neighborBlk = this->coveringLeafBlock( leafBlk.faceNeighbor(i) );
            if(neighborBlk != BlockIndex::invalid_index() )
            {
                InOutBlockData& neighborData = (*this)[neighborBlk];
                if (! neighborData.isColored() )
                {
                    QUEST_OCTREE_DEBUG_LOG_IF(
                            DEBUG_BLOCK_1 == neighborBlk || DEBUG_BLOCK_2 == neighborBlk
                             || DEBUG_BLOCK_1 == leafBlk || DEBUG_BLOCK_2 == leafBlk
                         , "Spreading color from block " << leafBlk << " with data " <<  leafData
                         << " and bounding box " << this->blockBoundingBox(leafBlk)
                         << " -- midpoint " << this->blockBoundingBox(leafBlk).centroid()
                         //
                         << "\n\t\t to " << neighborBlk << " with data " << neighborData
                         << " and bounding box " << this->blockBoundingBox(neighborBlk)
                         << " -- midpoint " << this->blockBoundingBox(neighborBlk).centroid()
                         );

                    switch(leafData.color())
                    {
                    case InOutBlockData::Black:
                        neighborData.setBlack();
                        break;
                    case InOutBlockData::White:
                        neighborData.setWhite();
                        break;
                    case InOutBlockData::Gray:
                        if( withinGrayBlock(this->blockBoundingBox(neighborBlk).centroid(),  leafBlk, leafData) )
                            neighborData.setBlack();
                        else
                            neighborData.setWhite();
                        break;
                    case InOutBlockData::Undetermined:
                        break;
                    }

                    QUEST_OCTREE_DEBUG_LOG_IF(neighborData.isColored()
                             && (DEBUG_BLOCK_1 == neighborBlk || DEBUG_BLOCK_2 == neighborBlk)
                          ,  "Neighbor block was colored -- " << neighborBlk
                             << " now has data " << neighborData);
                }
            }
        }
    }

    return isColored;
}


template<int DIM>
typename InOutOctree<DIM>::VertexIndex InOutOctree<DIM>::leafVertex(const BlockIndex& leafBlk, const InOutBlockData&  leafData) const
{
    if( m_generationState >= INOUTOCTREE_ELEMENTS_INSERTED)
    {
        SLIC_ASSERT( leafData.hasData() );
        return m_grayLeafToVertexRelationLevelMap[leafBlk.level()][leafData.dataIndex()][0];
    }
    else
    {
        return leafData.dataIndex();
    }
}

template<int DIM>
typename InOutOctree<DIM>::TriangleIndexSet InOutOctree<DIM>::leafTriangles(const BlockIndex& leafBlk, const InOutBlockData&  leafData) const
{
    SLIC_ASSERT( m_generationState >= INOUTOCTREE_ELEMENTS_INSERTED && leafData.hasData() );
    return m_grayLeafToElementRelationLevelMap[leafBlk.level()][leafData.dataIndex()];
}



template<int DIM>
bool InOutOctree<DIM>::withinGrayBlock(const SpacePt & pt, const BlockIndex& leafBlk, const InOutBlockData&  leafData) const
{
    SLIC_ASSERT( leafData.color() == InOutBlockData::Gray );
    SLIC_ASSERT( leafData.hasData() );

    VertexIndex vIdx = leafVertex(leafBlk, leafData);
    const SpaceVector vec( m_meshWrapper.vertexPosition( vIdx), pt );

    // If signs of all dot products are positive (negative) point is outside (inside)
    // Otherwise, short circuit for more complicated test below
    int sgnCount = 0;
    TriangleIndexSet triSet = leafTriangles(leafBlk, leafData);
    const int numTris = triSet.size();

    std::vector<SpaceTriangle> spaceTris;
    spaceTris.reserve( numTris);

    std::vector<SpaceVector>   unitNorms;
    unitNorms.reserve( numTris);

    QUEST_OCTREE_DEBUG_LOG_IF( DEBUG_BLOCK_1 == leafBlk || DEBUG_BLOCK_2 == leafBlk
          , "Within gray block " << leafBlk
          << "\n\t\t -- data " << leafData
          << "\n\t\t -- bounding box " << this->blockBoundingBox(leafBlk)
          << "\n\t\t -- testing point " << pt
          <<"\n\t\t -- block vertex is " << m_meshWrapper.vertexPosition( vIdx)
          << " -- vec is " << vec << " -- index " << vIdx
          <<"\n\t\t -- block triangles are " << triSet
        );

    for(int i=0; i< numTris; ++i)
    {
        TriangleIndex tIdx = triSet[i];
        spaceTris[i] = m_meshWrapper.trianglePositions(tIdx);
        unitNorms[i] = spaceTris[i].normal().unitVector();
        sgnCount += (vec.dot( unitNorms[i] )  >= 0 ) ? 1 : -1;
    }

    if(sgnCount == numTris)             // outside
        return false;
    else if(sgnCount == -numTris)       // inside
        return true;


    // Else, pt is within some half-spaces and outside others
    // Find the min distance to the closest triangle within leafBlk to pt
    // Note: we must account for the possibility that the closest point in on an edge or a vertex of a triangle.
    //       As a simple fix, we take the average of the normals.
    //       It might be more accurate to take the angle-weighted average as discussed in
    //          Baerentzen TVCG 11:3 (2005) -- `Signed distance computation using the angle weighted pseudonomral'

    std::vector<double> sqDists;
    sqDists.reserve( numTris);
    double minDistSq = std::numeric_limits<double>::max();

    for(int i=0; i< numTris; ++i)
    {
        sqDists[i] = squared_distance(pt, spaceTris[i] );
        if(sqDists[i] < minDistSq )
            minDistSq = sqDists[i];
    }

    SpaceVector norm;
    for(int i=0; i< numTris; ++i)
    {
        if( axom::utilities::isNearlyEqual(minDistSq, sqDists[i]))
            norm += unitNorms[i];
    }

    return vec.dot( norm.unitVector() ) >=0 ? false : true;
}


template<int DIM>
void InOutOctree<DIM>::updateSurfaceMeshVertices()
{
    // Create a map from old vertex indices to new vertex indices
    MeshVertexSet origVerts( m_meshWrapper.numMeshVertices() );
    VertexIndexMap vertexIndexMap( &origVerts, MeshWrapper::NO_VERTEX);

    // Generate unique indices for new mesh vertices
    int uniqueVertexCounter = 0;
    for(int i=0; i< origVerts.size(); ++i)
    {
        // Find the block and its indexed vertex in the octree
        BlockIndex leafBlock = this->findLeafBlock( m_meshWrapper.getMeshVertexPosition(i) );
        SLIC_ASSERT( (*this)[leafBlock].hasData() );
        VertexIndex vInd = (*this)[ leafBlock ].dataIndex();

        // If the indexed vertex doesn't have a new id, give it one
        if(vertexIndexMap[vInd] == MeshWrapper::NO_VERTEX)
            vertexIndexMap[vInd] = uniqueVertexCounter++;

        // If this is not the indexed vertex of the block, set the new index
        if(vInd != i)
            vertexIndexMap[i] = vertexIndexMap[vInd];
    }

    // Use the index map to reindex the mesh verts and elements
    m_meshWrapper.reindexMesh( uniqueVertexCounter, vertexIndexMap);

    // Update the octree leaf vertex IDs to the new mesh IDs
    // and create the map from the new vertices to their octree blocks
    m_vertexToBlockMap = VertexBlockMap(&m_meshWrapper.vertexSet());
    for(int i = 0; i< m_meshWrapper.numMeshVertices(); ++i)
    {
        const SpacePt& pos = m_meshWrapper.vertexPosition(i);
        BlockIndex leafBlock = this->findLeafBlock(pos);
        SLIC_ASSERT( this->isLeaf(leafBlock) && (*this)[leafBlock].hasData() );

        (*this)[ leafBlock ].setData(i);
        m_vertexToBlockMap[i] = leafBlock;
    }
}

template<int DIM>
bool InOutOctree<DIM>::allTrianglesIncidentInCommonVertex(const BlockIndex& leafBlock
                                                         , DynamicGrayBlockData& leafData) const
{
    bool shareCommonVert = false;

    VertexIndex commonVert = leafData.vertexIndex();

    const int numTris = leafData.numTriangles();
    const DynamicGrayBlockData::TriangleList& tris = leafData.triangles();

    if( blockIndexesVertex(commonVert, leafBlock) )
    {
        // This is a leaf node containing the indexed vertex
        // Loop through the triangles and check that all are incident with this vertex
        for(int i=0; i< numTris; ++i)
        {
            if( !m_meshWrapper.incidentInVertex(m_meshWrapper.triangleVertexIndices(tris[i]), commonVert) )
            {
                return false;
            }
        }
        shareCommonVert = true;
    }
    else
    {
        SLIC_ASSERT(numTris > 0);
        switch(numTris)
        {
        case 1:
            /// Choose an arbitrary vertex from this triangle
            commonVert = m_meshWrapper.triangleVertexIndices(tris[0])[0];
            shareCommonVert = true;
            break;
        case 2:
            /// Find a vertex that both triangles share
            shareCommonVert = m_meshWrapper.haveSharedVertex(tris[0],tris[1], commonVert);
            break;
        default :   // numTris >= 3
            /// Find a vertex that the first three triangles share
            shareCommonVert = m_meshWrapper.haveSharedVertex(tris[0],tris[1],tris[2], commonVert);

            /// Check that all other triangles have this vertex
            for(int i=3; shareCommonVert && i< numTris; ++i)
            {
                if( !m_meshWrapper.incidentInVertex( m_meshWrapper.triangleVertexIndices(tris[i]), commonVert) )
                        shareCommonVert = false;
            }
            break;
        }

        if(shareCommonVert)
            leafData.setVertex( commonVert );
    }

    return shareCommonVert;
}


template<int DIM>
bool InOutOctree<DIM>::within(const SpacePt& pt) const
{
    if( this->boundingBox().contains(pt) )
    {
        const BlockIndex block = this->findLeafBlock(pt);
        const InOutBlockData& data = (*this)[block];

        switch( data.color() )
        {
        case InOutBlockData::Black:
            return true;
        case InOutBlockData::White:
            return false;
        case InOutBlockData::Gray:
            return withinGrayBlock( pt, block, data);
        case InOutBlockData::Undetermined:
            SLIC_ASSERT_MSG(false, "Error -- All leaf blocks must have a color. "
                        << " The color of leafBlock " << block
                        << " was 'Undetermined' when querying point " << pt  );
            break;
        }
    }

    return false;
}


template<int DIM>
void InOutOctree<DIM>::printOctreeStats() const
{
    typedef typename OctreeBaseType::OctreeLevelType LeavesLevelMap;
    typedef typename OctreeBaseType::LevelMapCIterator LeavesIterator;
    typedef axom::slam::Map<int> LeafCountMap;

    LeafCountMap levelBlocks( &this->m_levels);
    LeafCountMap levelLeaves( &this->m_levels);
    LeafCountMap levelLeavesWithVert( &this->m_levels);
    LeafCountMap levelTriangleRefCount( &this->m_levels);

    LeafCountMap levelWhiteBlockCount( &this->m_levels);
    LeafCountMap levelBlackBlockCount( &this->m_levels);
    LeafCountMap levelGrayBlockCount( &this->m_levels);

    int totalBlocks = 0;
    int totalLeaves = 0;
    int totalLeavesWithVert = 0;
    int totalTriangleRefCount = 0;
    int totalWhiteBlocks = 0;
    int totalBlackBlocks = 0;
    int totalGrayBlocks = 0;

    // Iterate through blocks -- count the numbers of internal and leaf blocks
    for(int lev=0; lev< this->m_levels.size(); ++lev)
    {
        const LeavesLevelMap& levelLeafMap = this->getOctreeLevel( lev );
        levelBlocks[lev] = levelLeafMap.numBlocks();
        levelLeaves[lev] = levelLeafMap.numLeafBlocks();
        levelLeavesWithVert[lev] = 0;
        levelTriangleRefCount[lev] = 0;
        levelWhiteBlockCount[lev] = 0;
        levelBlackBlockCount[lev] = 0;
        levelGrayBlockCount[lev] = 0;

        for(LeavesIterator it=levelLeafMap.begin(), itEnd = levelLeafMap.end(); it != itEnd; ++it)
        {
            const InOutBlockData& blockData = *it;
            BlockIndex block(it.pt(), lev);

            if(blockData.isLeaf())
            {
                if(blockData.hasData())
                {
                    ++levelLeavesWithVert[ lev ];

                    if(m_generationState >= INOUTOCTREE_ELEMENTS_INSERTED)
                    {
                        levelTriangleRefCount[ lev ] += leafTriangles(block, blockData).size();
                    }
                }
                if(m_generationState >= INOUTOCTREE_LEAVES_COLORED)
                {
                    switch(blockData.color())
                    {
                    case InOutBlockData::Black:      ++levelBlackBlockCount[lev]; break;
                    case InOutBlockData::White:      ++levelWhiteBlockCount[lev]; break;
                    case InOutBlockData::Gray:       ++levelGrayBlockCount[lev];  break;
                    case InOutBlockData::Undetermined:                            break;
                    }
                }
            }
        }

        totalBlocks += levelBlocks[lev];
        totalLeaves += levelLeaves[lev];
        totalLeavesWithVert += levelLeavesWithVert[lev];
        totalTriangleRefCount += levelTriangleRefCount[ lev ];
        totalWhiteBlocks += levelWhiteBlockCount[lev];
        totalBlackBlocks += levelBlackBlockCount[lev];
        totalGrayBlocks  += levelGrayBlockCount[lev];
    }


    std::stringstream octreeStatsStr;
    octreeStatsStr << "*** " << (m_generationState >= INOUTOCTREE_ELEMENTS_INSERTED ? "PM" : "PR")
                   << " octree summary *** \n";

    for(int lev=0; lev< this->m_levels.size(); ++lev)
    {
        if(levelBlocks[lev] > 0)
        {
            int percentWithVert = (levelLeaves[lev] > 0)
                    ? (100. * levelLeavesWithVert[lev]) / levelLeaves[lev]
                    : 0;

            octreeStatsStr << "\t Level " << lev
                           << " has " << levelBlocks[lev] << " blocks -- "
                           <<  levelBlocks[lev] - levelLeaves[lev] << " internal; "
                           <<  levelLeaves[lev] << " leaves "
                           << " (" <<   percentWithVert  << "% w/ vert); ";
            if(m_generationState >= INOUTOCTREE_LEAVES_COLORED)
            {
                octreeStatsStr << " Leaves with colors -- B,W,G ==> " << levelBlackBlockCount[lev]
                               << "," << levelWhiteBlockCount[lev]
                               << "," << levelGrayBlockCount[lev]
                               << " and " << levelTriangleRefCount[lev] << " triangle references.";
            }
            //octreeStatsStr <<"Hash load factor: " << this->m_leavesLevelMap[ lev ].load_factor()
            //                                      << " -- max lf: " << this->m_leavesLevelMap[ lev ].max_load_factor();

            octreeStatsStr <<"\n";
        }
    }


    double meshNumTriangles = m_meshWrapper.numMeshElements();
    int percentWithVert = (100. * totalLeavesWithVert) / totalLeaves;
    octreeStatsStr<<"  Mesh has " << m_meshWrapper.numMeshVertices() << " vertices."
                 <<"\n  Octree has " << totalBlocks << " blocks; "
                 <<  totalBlocks - totalLeaves << " internal; "
                 <<  totalLeaves << " leaves "
                 << " (" <<   percentWithVert  << "% w/ vert); ";
    if(m_generationState >= INOUTOCTREE_ELEMENTS_INSERTED)
    {
        octreeStatsStr<<" \n\t There were " << totalTriangleRefCount << " triangle references "
                 <<" (avg. " << ( totalTriangleRefCount / meshNumTriangles ) << " refs per triangle).";
    }
    if(m_generationState >= INOUTOCTREE_LEAVES_COLORED)
    {
        octreeStatsStr<<"\n Colors B,W,G ==> " << totalBlackBlocks
                 << "," << totalWhiteBlocks
                 << "," << totalGrayBlocks
                 <<"\n";
    }

    SLIC_INFO( octreeStatsStr.str() );

    if(m_generationState >= INOUTOCTREE_ELEMENTS_INSERTED)
    {
        if(DEBUG_TRI_IDX >= 0)
        {
             dumpMeshVTK("triangle", DEBUG_TRI_IDX, this->root(), this->blockBoundingBox(this->root()), true);


             TriVertIndices tv = m_meshWrapper.triangleVertexIndices(DEBUG_TRI_IDX);
             for(int i=0; i< 3; ++i)
             {
                 BlockIndex blk = m_vertexToBlockMap[ tv[i] ];
                 dumpMeshVTK("triangleVertexBlock", DEBUG_TRI_IDX, blk, this->blockBoundingBox(blk), false);
             }

        }

        typedef axom::slam::Map<int> TriCountMap;
        typedef axom::slam::Map<int> CardinalityVTMap;

        TriCountMap triCount( &m_meshWrapper.elementSet());
        for(int lev=0; lev< this->m_levels.size(); ++lev)
        {
            const LeavesLevelMap& levelLeafMap = this->getOctreeLevel( lev );
            levelBlocks[lev] = levelLeafMap.numBlocks();
            for(LeavesIterator it=levelLeafMap.begin(), itEnd = levelLeafMap.end(); it != itEnd; ++it)
            {
                const InOutBlockData& leafData = *it;
                if(leafData.isLeaf() && leafData.hasData())
                {
                    BlockIndex blk(it.pt(), lev);
                    TriangleIndexSet tris = leafTriangles(blk, leafData);
                    for(int i = 0; i < tris.size(); ++i)
                    {
                        ++triCount[ tris[i]];

                        if(DEBUG_TRI_IDX == tris[i] )
                        {
                            dumpMeshVTK("triangleIntersect", tris[i], blk, this->blockBoundingBox(blk), false);
                        }

                    }
                }
            }
        }


        // Generate and output a histogram of the bucket counts on a lg-scale
        typedef std::map<int,int> LogHistogram;
        LogHistogram triCountHist;        // Create histogram of edge lengths (log scale)

        typedef axom::primal::BoundingBox<double,1> MinMaxRange;
        typedef MinMaxRange::PointType LengthType;

        typedef std::map<int,MinMaxRange> LogRangeMap;
        LogRangeMap triCountRange;

        for ( int i=0; i < m_meshWrapper.numMeshElements(); ++i )
        {
            LengthType count( triCount[i] );
            int expBase2;
            std::frexp( triCount[i], &expBase2);
            triCountHist[ expBase2 ]++;
            triCountRange[ expBase2 ].addPoint( count );
        }

        std::stringstream triCountStr;
        triCountStr<<"\tTriangle index count (lg-arithmic bins for number of references per triangle):";
        for(LogHistogram::const_iterator it = triCountHist.begin()
                ; it != triCountHist.end()
                ; ++it)
        {
            triCountStr << "\n\t exp: " << it->first
                        <<"\t count: " << (it->second)
                        <<"\tRange: " << triCountRange[it->first];
        }
        SLIC_INFO(triCountStr.str());

        // Generate and output histogram of VT relation
        CardinalityVTMap cardVT(&m_meshWrapper.vertexSet());
        for ( int i=0; i < m_meshWrapper.numMeshElements(); ++i )
        {
            TriVertIndices tvRel = m_meshWrapper.triangleVertexIndices(i);
            cardVT[tvRel[0]]++;
            cardVT[tvRel[1]]++;
            cardVT[tvRel[2]]++;
        }

        typedef std::map<int,int> LinHistogram;
        LinHistogram vtCardHist;
        for ( int i=0; i < m_meshWrapper.numMeshVertices(); ++i )
        {
            LengthType count( cardVT[i] );
            vtCardHist[ cardVT[i] ]++;
        }

        std::stringstream vtCartStr;
        vtCartStr<<"\tCardinality VT relation histogram (linear): ";
        for(LinHistogram::const_iterator it = vtCardHist.begin()
                ; it != vtCardHist.end()
                ; ++it)
        {
            vtCartStr << "\n\t exp: " << it->first
                        <<"\t count: " << (it->second);
        }
        SLIC_INFO(vtCartStr.str());

//        // Add field to the triangle mesh
//        mint::FieldData* CD = m_surfaceMesh->getCellFieldData();
//        CD->addField( new mint::FieldVariable< TriangleIndex >("blockCount", meshTris.size()) );
//
//        int* blockCount = CD->getField( "blockCount" )->getIntPtr();
//
//        SLIC_ASSERT( blockCount != AXOM_NULLPTR );
//
//        for ( int i=0; i < meshTris.size(); ++i ) {
//            blockCount[i] = triCount[i];
//        }
//
//        // Add field to the triangle mesh
//        mint::FieldData* ND = m_surfaceMesh->getNodeFieldData();
//        ND->addField( new mint::FieldVariable< int >("vtCount", m_vertexSet.size()) );
//
//        int* vtCount = ND->getField( "vtCount" )->getIntPtr();
//
//        SLIC_ASSERT( vtCount != AXOM_NULLPTR );
//
//        for ( int i=0; i < m_vertexSet.size(); ++i ) {
//            vtCount[i] = cardVT[i];
//        }
    }

}



template<int DIM>
void InOutOctree<DIM>::checkAllLeavesColoredAtLevel(int AXOM_DEBUG_PARAM(level)) const
{
#ifdef AXOM_DEBUG
    typedef typename OctreeBaseType::OctreeLevelType LeavesLevelMap;
    typedef typename OctreeBaseType::LevelMapCIterator LeavesIterator;

  const LeavesLevelMap& levelLeafMap = this->getOctreeLevel( level );
  for(LeavesIterator it=levelLeafMap.begin(), itEnd = levelLeafMap.end(); it != itEnd; ++it)
  {
      if( ! it->isLeaf() )
          continue;

      SLIC_ASSERT_MSG( it->isColored()
                      , "Error after coloring level " << level
                      << " leaf block " << BlockIndex(it.pt(), level)
                      << " was not colored."
      );
  }
 #endif
}

template<int DIM>
void InOutOctree<DIM>::checkValid() const
{
#ifdef AXOM_DEBUG
    typedef typename OctreeBaseType::OctreeLevelType LeavesLevelMap;
    typedef typename OctreeBaseType::LevelMapCIterator LeavesIterator;

    SLIC_ASSERT_MSG( m_generationState >= INOUTOCTREE_VERTICES_INSERTED
                   , "InOutOctree::checkValid assumes that we have already inserted "
                     << "the vertices into the octree.");

    // Iterate through the tree
    // Internal blocks should not have associated vertex data
    // Leaf block consistency depends on 'color'
    //      Black or White blocks should not have vertex data
    //      Gray blocks should have a vertex reference; it may or may not be located within the block
    //          The sum of vertices located within a block should equal the number of mesh vertices.
    //      Gray blocks should have one or more triangles.
    //          All triangles should be incident in a common vertex -- which equals the indexed vertex, if it exists.

    std::string treeTypeStr = (m_generationState >= INOUTOCTREE_ELEMENTS_INSERTED)
                ? "PM" : "PR";

    SLIC_DEBUG("Inside InOutOctree::checkValid() to verify state of " << treeTypeStr << " octree.");


    if(m_generationState > INOUTOCTREE_VERTICES_INSERTED )
    {
      SLIC_DEBUG("--Checking that each vertex is in a leaf block of the tree.");
      const int numVertices = m_meshWrapper.numMeshVertices();
      for(int i=0; i< numVertices; ++i)
      {
        const SpacePt& pos = m_meshWrapper.vertexPosition(i);
        BlockIndex vertBlock = this->findLeafBlock(pos);
        const InOutBlockData& leafData = (*this)[vertBlock];
        VertexIndex vertInBlock = leafVertex(vertBlock, leafData);

        // Check that we can find the leaf block indexing each vertex from its position
        SLIC_ASSERT_MSG( leafData.hasData() &&  vertInBlock == i
                     ,  " Vertex " << i << " at position " << pos
                     << " \n\t was not indexed in block " << vertBlock
                     << " with bounding box " << this->blockBoundingBox(vertBlock)
                     << " ( point is" << (this->blockBoundingBox(vertBlock).contains(pos)? "" : " NOT" )
                     << " contained in block )."
                     << " \n\n *** \n Leaf data: " << leafData
                     << " \n ***"
        );

        // Check that our cached value of the vertex's block is accurate
        SLIC_ASSERT_MSG( vertBlock == m_vertexToBlockMap[i]
                   , "Cached block for vertex " << i << " differs from its found block. "
                   << "\n\t -- cached block "<< m_vertexToBlockMap[i]
                   << "-- is leaf? " << ( (*this)[ m_vertexToBlockMap[i]].isLeaf() ? "yes" : "no")
                   << "\n\t -- actual block " << vertBlock
                   << "-- is leaf? " << ( (*this)[ vertBlock].isLeaf() ? "yes" : "no")
                   << "\n\t -- vertex set size: " << numVertices
                    );
      }
    }

    if(m_generationState >= INOUTOCTREE_ELEMENTS_INSERTED)
    {
        SLIC_DEBUG("--Checking that each triangle is referenced by the leaf blocks containing its vertices.");

        const int numTriangles = m_meshWrapper.numMeshElements();
        for(int tIdx=0; tIdx< numTriangles; ++tIdx)
        {
            TriVertIndices tvRel = m_meshWrapper.triangleVertexIndices( tIdx );
            for(int j=0; j< tvRel.size();++j)
            {
                VertexIndex vIdx = tvRel[j];
                BlockIndex vertBlock = m_vertexToBlockMap[vIdx];
                const InOutBlockData& leafData = (*this)[vertBlock];

                // Check that this triangle is referenced here.
                bool foundTriangle = false;
                TriangleIndexSet leafTris = leafTriangles(vertBlock, leafData);
                for(int k=0; !foundTriangle && k< leafTris.size(); ++k)
                {
                    if( leafTris[k] == tIdx)
                        foundTriangle = true;
                }

                SLIC_ASSERT_MSG( foundTriangle
                               , "Did not find triangle with index " << tIdx << " and vertices" << tvRel
                               << " in block " << vertBlock << " containing vertex " << vIdx
                               << " \n\n *** \n Leaf data: " << leafData
                               << " \n\t Triangles in block? " << leafTris
                               << " \n ***"
                               );
            }
        }

        // Check that internal blocks have no triangle / vertex
        //       and leaf blocks satisfy the conditions above
        SLIC_DEBUG("--Checking that internal blocks have no data, and that leaves satisfy all PM conditions");
        for(int lev=0; lev< this->m_levels.size(); ++lev)
        {
            const LeavesLevelMap& levelLeafMap = this->getOctreeLevel(lev);
            for(LeavesIterator it=levelLeafMap.begin(), itEnd = levelLeafMap.end(); it != itEnd; ++it)
            {
                const BlockIndex block(it.pt(), lev);
                const InOutBlockData& data = *it;

                if( !data.isLeaf() )
                {
                    SLIC_ASSERT( !data.hasData() );
                }
                else // leaf block
                {
                    if( data.hasData())
                    {
                        VertexIndex vIdx = leafVertex(block, data);
                        TriangleIndexSet triSet = leafTriangles(block,data);
                        for( int i  = 0; i< triSet.size(); ++i)
                        {
                            TriangleIndex tIdx = triSet[i];

                            // Check that vIdx is one of this triangle's vertices
                            TriVertIndices tvRel = m_meshWrapper.triangleVertexIndices( tIdx );
                            SLIC_ASSERT_MSG( m_meshWrapper.incidentInVertex(tvRel, vIdx)
                                             , "All triangles in a gray block must be incident on a common vertex,"
                                             << " but triangles " << tIdx << " with vertices " << tvRel
                                             << " in block " << block << " is not incident in vertex " << vIdx);

                            // Check that this triangle intersects the bounding box of the block
                            GeometricBoundingBox blockBB = this->blockBoundingBox(block);
                            SLIC_ASSERT_MSG( blockIndexesElementVertex(tIdx, block)
                                             || intersect( m_meshWrapper.trianglePositions(tIdx), blockBB)
                                           , "Triangle " << tIdx << " was indexed in block " << block
                                            << " but it does not intersect the block."
                                            << "\n\tBlock bounding box: " << blockBB
                                            << "\n\tTriangle positions: " << m_meshWrapper.trianglePositions(tIdx)
                                            << "\n\tTriangle vertex indices: " << tvRel
                                            << "\n\tLeaf vertex is: " << vIdx
                                            << "\n\tLeaf triangles: " << triSet
                                            << "(" << triSet.size() <<")"
                            );
                        }
                    }
                }

            }
        }
    }

    if(m_generationState >= INOUTOCTREE_LEAVES_COLORED)
    {
        SLIC_DEBUG("--Checking that all leaves have a color -- black, white and gray");
        for(int lev=0; lev< this->m_levels.size(); ++lev)
        {
            checkAllLeavesColoredAtLevel(lev);
        }

        SLIC_DEBUG("--Checking that inside blocks do not neighbor outside blocks");
        for(int lev=this->maxLeafLevel()-1; lev >= 0; --lev)
        {
            const LeavesLevelMap& levelLeafMap = this->getOctreeLevel( lev );
            for(LeavesIterator it=levelLeafMap.begin(), itEnd = levelLeafMap.end(); it != itEnd; ++it)
            {
                const BlockIndex block(it.pt(), lev);
                const InOutBlockData& data = *it;

                if(data.isLeaf() && data.color() != InOutBlockData::Gray)
                {
                    // Leaf does not yet have a color... try to find its color from same-level face neighbors
                    for(int i=0; i< block.numFaceNeighbors(); ++i)
                    {
                        BlockIndex neighborBlk = this->coveringLeafBlock( block.faceNeighbor(i) );
                        if(neighborBlk != BlockIndex::invalid_index() )
                        {
                            const InOutBlockData& neighborData = (*this)[neighborBlk];
                            switch(neighborData.color())
                            {
                            case InOutBlockData::Black:             // intentional fallthrough
                            case InOutBlockData::White:
                                SLIC_CHECK_MSG(data.color() == neighborData.color()
                                                , "Problem at block " << block << " with data " << data
                                                <<" --- neighbor " << neighborBlk << " with data " << neighborData
                                                <<". Neighboring leaves that are not gray must have the same color."
                                                );
                                break;
                            case InOutBlockData::Gray:              // intentional fallthrough
                            case InOutBlockData::Undetermined:
                                break;
                            }
                        }
                    }



                }
            }
        }

    }

#endif
    SLIC_DEBUG("done.");

}
template<int DIM>
void InOutOctree<DIM>::dumpOctreeMeshVTK( const std::string& name) const
{
  #ifdef DUMP_VTK_MESH
    typedef typename OctreeBaseType::OctreeLevelType LeavesLevelMap;
    typedef typename OctreeBaseType::LevelMapCIterator LeavesIterator;


    DebugMesh* debugMesh= new DebugMesh(3);
    std::stringstream fNameStr;

    int totalLeaves = 0;

    const bool hasTriangles = (m_generationState >= INOUTOCTREE_ELEMENTS_INSERTED);
    const bool hasColors = (m_generationState >= INOUTOCTREE_LEAVES_COLORED);

    // Iterate through blocks -- count the numbers of leaves
    for(int lev=0; lev< this->m_levels.size(); ++lev)
    {
        const LeavesLevelMap& levelLeafMap = this->getOctreeLevel( lev );
        for(LeavesIterator it=levelLeafMap.begin(), itEnd = levelLeafMap.end(); it != itEnd; ++it)
        {
            if( it->isLeaf())
                totalLeaves++;
        }
    }
    SLIC_INFO("Dump vtk:: Octree has " << totalLeaves << " leaves.");

    // Create a small lookup table to map block colors to ints
    std::map< InOutBlockData::LeafColor, int> colorMap;
    colorMap[ InOutBlockData::White ] = -1;
    colorMap[ InOutBlockData::Gray ] = 0;
    colorMap[ InOutBlockData::Black ] = 1;

    // Allocate Slam Maps for the field data
    axom::slam::PositionSet leafSet(totalLeaves);
    typedef axom::slam::Map<VertexIndex> LeafVertMap;
    typedef axom::slam::Map<int> LeafIntMap;
    LeafVertMap leafVertID(&leafSet);
    LeafVertMap leafVertID_unique(&leafSet);
    LeafIntMap leafTriCount(&leafSet);
    LeafIntMap leafColors(&leafSet);
    LeafIntMap leafLevel(&leafSet);


    // Iterate through blocks -- and set the field data
    int leafCount = 0;
    for(int lev=0; lev< this->m_levels.size(); ++lev)
    {
        const LeavesLevelMap& levelLeafMap = this->getOctreeLevel( lev );
        for(LeavesIterator it=levelLeafMap.begin(), itEnd = levelLeafMap.end(); it != itEnd; ++it)
        {
            if( ! it->isLeaf())
                continue;

            // Add a hex
            BlockIndex block( it.pt(), lev);
            const InOutBlockData& leafData = *it;
            GeometricBoundingBox blockBB = this->blockBoundingBox(block);

            int vStart = debugMesh->getMeshNumberOfNodes();
            debugMesh->insertNode( blockBB.getMin()[0], blockBB.getMin()[1], blockBB.getMin()[2]);
            debugMesh->insertNode( blockBB.getMax()[0], blockBB.getMin()[1], blockBB.getMin()[2]);
            debugMesh->insertNode( blockBB.getMax()[0], blockBB.getMax()[1], blockBB.getMin()[2]);
            debugMesh->insertNode( blockBB.getMin()[0], blockBB.getMax()[1], blockBB.getMin()[2]);

            debugMesh->insertNode( blockBB.getMin()[0], blockBB.getMin()[1], blockBB.getMax()[2]);
            debugMesh->insertNode( blockBB.getMax()[0], blockBB.getMin()[1], blockBB.getMax()[2]);
            debugMesh->insertNode( blockBB.getMax()[0], blockBB.getMax()[1], blockBB.getMax()[2]);
            debugMesh->insertNode( blockBB.getMin()[0], blockBB.getMax()[1], blockBB.getMax()[2]);

            int data[8];
            for(int i=0; i< 8; ++i)
                data[i] = vStart + i;

            debugMesh->insertCell( data, mint::LINEAR_HEX, 8);

            int vIdx = leafData.hasData() ? leafVertex(block, leafData) :  MeshWrapper::NO_VERTEX;
            leafVertID[leafCount] = vIdx;

            leafLevel[leafCount] = lev;

            if(hasTriangles)
            {
                leafVertID_unique[leafCount] = blockIndexesVertex(vIdx, block) ? vIdx : MeshWrapper::NO_VERTEX;
                leafTriCount[leafCount] = leafData.hasData()? leafTriangles(block,leafData).size() : 0;
            }

            if(hasColors)
            {
                leafColors[leafCount] = colorMap[ leafData.color()];
            }

            leafCount++;
        }
    }

    if(totalLeaves > 0)
    {
        // Add the fields to the mesh
        mint::FieldData* CD = debugMesh->getCellFieldData();

        CD->addField( new mint::FieldVariable< VertexIndex >("vertID", leafSet.size()) );
        VertexIndex* vertID = CD->getField( "vertID" )->getIntPtr();
        SLIC_ASSERT( vertID != AXOM_NULLPTR );

        CD->addField( new mint::FieldVariable< VertexIndex >("level", leafSet.size()) );
        VertexIndex* lLevel = CD->getField( "level" )->getIntPtr();
        SLIC_ASSERT( lLevel != AXOM_NULLPTR );

        for ( int i=0; i < leafSet.size(); ++i )
        {
            vertID[i] = leafVertID[i];
            lLevel[i] = leafLevel[i];
        }

        if(hasTriangles)
        {
            CD->addField( new mint::FieldVariable< VertexIndex >("uniqVertID", leafSet.size()) );
            VertexIndex* uniqVertID = hasTriangles ? CD->getField( "uniqVertID" )->getIntPtr(): AXOM_NULLPTR;
            SLIC_ASSERT( uniqVertID != AXOM_NULLPTR );

            CD->addField( new mint::FieldVariable< int >("triCount", leafSet.size()) );
            int* triCount = hasTriangles ? CD->getField( "triCount" )->getIntPtr() : AXOM_NULLPTR;
            SLIC_ASSERT( triCount != AXOM_NULLPTR );

            for ( int i=0; i < leafSet.size(); ++i )
            {
                uniqVertID[i] = leafVertID_unique[i];
                triCount[i] = leafTriCount[i];
            }
        }

        if(hasColors)
        {
            CD->addField( new mint::FieldVariable< int >("colors", leafSet.size()) );
            int* colors = hasColors ? CD->getField( "colors" )->getIntPtr() : AXOM_NULLPTR;
            SLIC_ASSERT( !hasColors || colors != AXOM_NULLPTR );
            for ( int i=0; i < leafSet.size(); ++i )
                colors[i] = leafColors[i];
        }

        debugMesh->toVtkFile(name);
    }

    delete debugMesh;
    debugMesh = AXOM_NULLPTR;
  #else
    // Do something with the parameters to avoid a warning about unused parameters
    AXOM_DEBUG_VAR(name);
  #endif
}


template<int DIM>
void InOutOctree<DIM>::dumpDifferentColoredNeighborsMeshVTK( const std::string& fName) const
{
  #ifdef DUMP_VTK_MESH
    typedef typename OctreeBaseType::OctreeLevelType LeavesLevelMap;
    typedef typename OctreeBaseType::LevelMapCIterator LeavesIterator;

    if( m_generationState <= INOUTOCTREE_LEAVES_COLORED)
    {
        SLIC_INFO("Need to generate octree colors before visualizing them.");
        return;
    }

    DebugMesh* debugMesh= new DebugMesh(3);

    std::vector<int> colorsField;
    std::vector<int> debugIdxField;

    std::map< InOutBlockData::LeafColor, int> colorMap;
    colorMap[ InOutBlockData::White ] = -1;
    colorMap[ InOutBlockData::Gray ] = 0;
    colorMap[ InOutBlockData::Black ] = 1;

#if defined(AXOM_USE_CXX11)
  typedef std::unordered_map<GridPt, int, PointHash<int> > GridIntMap;
#else
  typedef boost::unordered_map<GridPt, int, PointHash<int> > GridIntMap;
#endif
  typedef typename GridIntMap::iterator GridIntIter;

  typedef axom::slam::Map<GridIntMap> LevelGridIntMap;
  LevelGridIntMap diffBlocks( &(this->m_levels) );

  // Iterate through the octree leaves looking for neigbor blocks with different labelings
  for(int lev=this->maxLeafLevel()-1; lev >= 0; --lev)
  {
    const LeavesLevelMap& levelLeafMap = this->getOctreeLevel( lev );
    for(LeavesIterator it=levelLeafMap.begin(), itEnd = levelLeafMap.end(); it != itEnd; ++it)
    {
      const BlockIndex block(it.pt(), lev);
      const InOutBlockData& data = *it;
      if(data.isLeaf() && data.color() != InOutBlockData::Gray)
      {
        bool addNeigbors = false;
        for(int i=0; i< block.numFaceNeighbors(); ++i)
        {
          BlockIndex neighborBlk = this->coveringLeafBlock( block.faceNeighbor(i) );
          if(neighborBlk != BlockIndex::invalid_index() )
          {
            const InOutBlockData& neighborData = (*this)[neighborBlk];
            switch(neighborData.color())
            {
            case InOutBlockData::Black:
            case InOutBlockData::White:
                if(data.color() != neighborData.color())
                {
                    diffBlocks[lev][block.pt()] = 1;
                    addNeigbors = true;
                }
                break;
            case InOutBlockData::Gray:
            case InOutBlockData::Undetermined:
                break;
            }
          }
        }
        if(addNeigbors)
        {
          for(int i=0; i< block.numFaceNeighbors(); ++i)
          {
            BlockIndex neighborBlk = this->coveringLeafBlock( block.faceNeighbor(i) );
            if(neighborBlk != BlockIndex::invalid_index() )
            {
                GridIntMap& levelBlocks = diffBlocks[ neighborBlk.level()];
                if( levelBlocks.find( neighborBlk.pt() ) == levelBlocks.end() )
                    levelBlocks[neighborBlk.pt()] = 0;
            }
          }
        }
      }
    }
  }

  // Add an offset to the debug blocks for easier identification
  if(this->inBounds(DEBUG_BLOCK_1) && this->hasBlock(DEBUG_BLOCK_1))
    diffBlocks[DEBUG_BLOCK_1.level()][ DEBUG_BLOCK_1.pt() ] += 10;
  if(this->inBounds(DEBUG_BLOCK_2) && this->hasBlock(DEBUG_BLOCK_2))
    diffBlocks[DEBUG_BLOCK_2.level()][ DEBUG_BLOCK_2.pt() ] += 20;

  for(int lev=this->maxLeafLevel()-1; lev >= 0; --lev)
  {
    GridIntMap& levelBlocks = diffBlocks[ lev ];
    for(GridIntIter it=levelBlocks.begin(), itEnd = levelBlocks.end(); it != itEnd; ++it)
    {
      const BlockIndex blk(it->first, lev);
      const InOutBlockData& blkData = (*this)[blk];

      GeometricBoundingBox blockBB = this->blockBoundingBox(blk);
      debugMesh->insertNode( blockBB.getMin()[0], blockBB.getMin()[1], blockBB.getMin()[2]);
      debugMesh->insertNode( blockBB.getMax()[0], blockBB.getMin()[1], blockBB.getMin()[2]);
      debugMesh->insertNode( blockBB.getMax()[0], blockBB.getMax()[1], blockBB.getMin()[2]);
      debugMesh->insertNode( blockBB.getMin()[0], blockBB.getMax()[1], blockBB.getMin()[2]);

      debugMesh->insertNode( blockBB.getMin()[0], blockBB.getMin()[1], blockBB.getMax()[2]);
      debugMesh->insertNode( blockBB.getMax()[0], blockBB.getMin()[1], blockBB.getMax()[2]);
      debugMesh->insertNode( blockBB.getMax()[0], blockBB.getMax()[1], blockBB.getMax()[2]);
      debugMesh->insertNode( blockBB.getMin()[0], blockBB.getMax()[1], blockBB.getMax()[2]);

      // Write the blocks with problematic neighbors
      int data[8];
      for(int i=0; i< 8; ++i)
        data[i] = colorsField.size()*8+i;
      debugMesh->insertCell( data, mint::LINEAR_HEX, 8);

      debugIdxField.push_back( it->second );
      colorsField.push_back( colorMap[ blkData.color() ]);
    }
  }

  if(!colorsField.empty())
  {
    mint::FieldData* CD = debugMesh->getCellFieldData();
    CD->addField( new mint::FieldVariable< int >("colors", colorsField.size()) );
    CD->addField( new mint::FieldVariable< int >("debugIdx", debugIdxField.size()) );
    int* colors = CD->getField( "colors" )->getIntPtr();
    int* debIdx = CD->getField( "debugIdx" )->getIntPtr();
    SLIC_ASSERT(colors != AXOM_NULLPTR );
    SLIC_ASSERT(debIdx != AXOM_NULLPTR );
    for ( std::size_t i=0; i < colorsField.size(); ++i )
    {
      debIdx[i] = debugIdxField[i];
      colors[i] = colorsField[i];
    }

    debugMesh->toVtkFile(fName);
  }

  delete debugMesh;
  debugMesh = AXOM_NULLPTR;
 #else
  // Do something with the parameters to avoid a warning about unused parameters
  AXOM_DEBUG_VAR(fName);
 #endif
}

template<int DIM>
void InOutOctree<DIM>::dumpMeshVTK( const std::string& name
                , int idx
                , const BlockIndex& block
                , const GeometricBoundingBox& blockBB
                , bool isTri
                ) const
{
  #ifdef DUMP_VTK_MESH
    DebugMesh* debugMesh= new DebugMesh(3);
    std::stringstream fNameStr;

    if(isTri)
    {
        SpaceTriangle triPos = m_meshWrapper.trianglePositions(idx);

        debugMesh->insertNode( triPos[0][0], triPos[0][1], triPos[0][2]);
        debugMesh->insertNode( triPos[1][0], triPos[1][1], triPos[1][2]);
        debugMesh->insertNode( triPos[2][0], triPos[2][1], triPos[2][2]);

        int verts[3] = {0,1,2};
        debugMesh->insertCell(verts, mint::LINEAR_TRIANGLE, 3);

        fNameStr << name << idx << ".vtk";

        SLIC_INFO("// Triangle " << idx );
        SLIC_INFO("TriangleType tri(PointType::make_point" << triPos[0]
                              << ", PointType::make_point" << triPos[1]
                              << ", PointType::make_point" << triPos[2]<<");" );

    }
    else
    {
        static std::map<std::string, int> counter;

        debugMesh->insertNode( blockBB.getMin()[0], blockBB.getMin()[1], blockBB.getMin()[2]);
        debugMesh->insertNode( blockBB.getMax()[0], blockBB.getMin()[1], blockBB.getMin()[2]);
        debugMesh->insertNode( blockBB.getMax()[0], blockBB.getMax()[1], blockBB.getMin()[2]);
        debugMesh->insertNode( blockBB.getMin()[0], blockBB.getMax()[1], blockBB.getMin()[2]);

        debugMesh->insertNode( blockBB.getMin()[0], blockBB.getMin()[1], blockBB.getMax()[2]);
        debugMesh->insertNode( blockBB.getMax()[0], blockBB.getMin()[1], blockBB.getMax()[2]);
        debugMesh->insertNode( blockBB.getMax()[0], blockBB.getMax()[1], blockBB.getMax()[2]);
        debugMesh->insertNode( blockBB.getMin()[0], blockBB.getMax()[1], blockBB.getMax()[2]);

        int data[8];
        for(int i=0; i< 8; ++i)
            data[i] = i;

        debugMesh->insertCell( data, mint::LINEAR_HEX, 8);

        fNameStr << name << idx << "_" << counter[name] << ".vtk";

        SLIC_INFO("// Block index " << block);
        SLIC_INFO("BoundingBoxType box"<<counter[name]<<"(PointType::make_point" << blockBB.getMin() << ", PointType::make_point" << blockBB.getMax()<<");" );
        counter[name]++;
    }

    debugMesh->toVtkFile(fNameStr.str());

    delete debugMesh;
    debugMesh = AXOM_NULLPTR;
  #else
    // Do something with the parameters to avoid a warning about unused parameters
    AXOM_DEBUG_VAR(name);
    AXOM_DEBUG_VAR(idx);
    AXOM_DEBUG_VAR(block);
    AXOM_DEBUG_VAR(blockBB);
    AXOM_DEBUG_VAR(isTri);
  #endif
}


} // end namespace quest 
} // end namespace axom 

#endif  // SPATIAL_OCTREE__HXX_
