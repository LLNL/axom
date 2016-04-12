#ifndef INOUT_OCTREE__HXX_
#define INOUT_OCTREE__HXX_

#include "common/ATKMacros.hpp"
#include "common/Timer.hpp"

#include "quest/BoundingBox.hpp"
#include "quest/Intersection.hpp"
#include "quest/Orientation.hpp"
#include "quest/Point.hpp"
#include "quest/Vector.hpp"
#include "quest/SquaredDistance.hpp"
#include "quest/fuzzy_compare.hpp"


#include "slic/slic.hpp"

#include "slam/Map.hpp"
#include "slam/RangeSet.hpp"
#include "slam/StaticConstantRelation.hpp"
#include "slam/StaticVariableRelation.hpp"


#include "quest/Mesh.hpp"
#include "quest/SpatialOctree.hpp"
#include "quest/Triangle.hpp"
#include "quest/Mesh.hpp"
#include "quest/UnstructuredMesh.hpp"
#include "quest/FieldData.hpp"
#include "quest/FieldVariable.hpp"


#include <vector> // For InOutLeafData triangle lists -- TODO replace with SLAM DynamicVariableRelation...
#include <ostream>
#include <set>
#include <iterator> // For back_inserter



#define DEBUG_VERT_IDX  -2
#define DEBUG_TRI_IDX   -2

#ifndef DUMP_VTK_MESH
//    #define DUMP_VTK_MESH
#endif

#ifndef DUMP_OCTREE_INFO
    #define DUMP_OCTREE_INFO 1
#endif


/**
 * \file
 *
 * \brief Defines an InOutOctree for containment queries on a surface.
 */


namespace quest
{

    class InOutBlockData
    {
        // Some internal constants for keeping tracking of the associated block
        // A block is a leaf block when its m_idx is not INTERNAL_BLOCK
        // Leaf blocks can be uncolored or colored (without additional data)
        //      or m_idx be the index of the data associated with a gray block
        enum { INTERNAL_BLOCK       = -1
             , LEAF_BLOCK_UNCOLORED = -2
             , LEAF_BLOCK_WHITE     = -3
             , LEAF_BLOCK_BLACK     = -4
        };

    public:
        enum LeafColor { Undetermined = -2, White = -1, Gray = 0, Black = 1};

    public:
        /**
         * \brief Default constructor for an InOutBlockData
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

        bool isLeaf() const { return m_idx != INTERNAL_BLOCK; }

        void setInternal()  { m_idx = INTERNAL_BLOCK; }

    public:     // Other functions

        void clear()
        {
            // No-op for now -- eventually, will need to do something about the index
        }

        bool hasData() const { return m_idx >= 0; }

        const int& dataIndex() const
        {
            //SLIC_ASSERT(hasData());
            return m_idx;
        }

        /**
         * \brief Sets the block as gray, and provides the index of its associated data
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

        void setBlack()
        {
            SLIC_ASSERT( isLeaf() );
            m_idx = LEAF_BLOCK_BLACK;
        }

        void setWhite()
        {
            SLIC_ASSERT( isLeaf() );
            m_idx = LEAF_BLOCK_WHITE;
        }

        void setData(int idx)
        {
            m_idx = idx;
        }

        void setUncoloredLeaf()
        {
            SLIC_ASSERT( isLeaf() );
            m_idx = LEAF_BLOCK_UNCOLORED;
        }

        /**
         * \brief Find the 'color' of this LeafBlock
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

        bool isColored() const
        {
            return color() != Undetermined;
        }

    private:
        int m_idx;
    };

    std::ostream& operator<<(std::ostream& os, const InOutBlockData& iob)
    {
        os << "InOutBlockData{"
           << "isLeaf: " << (iob.isLeaf() ? "yes" : "no");


        os<<", dataIndex: ";
        if(iob.hasData())
            os << "<no data>";
        else
            os << iob.dataIndex();

        os << "}";

        return os;
    }

    /**
     * \brief DataType for InOut octree leaf nodes
     * Contains a reference to a vertex ID in the mesh
     * (and has a unique ID from its parent class)
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
        bool isLeaf() const { return m_isLeaf; }
        void setLeafFlag(bool isLeaf) { m_isLeaf = isLeaf; }

    public: // Functions related to the associated vertex

        /**
         * \brief Checks whether there is a vertex associated with this leaf
         */
        bool hasVertex() const { return m_vertIndex >= 0; }

        /**
         * \brief Sets the vertex associated with this leaf
         */
        void setVertex(VertexIndex vInd) { m_vertIndex = vInd; }

        void clearVertex() { m_vertIndex = NO_VERTEX; }

        /**
         * \brief Accessor for the index of the vertex associated with this leaf
         */
        VertexIndex& vertexIndex() { return m_vertIndex; }

        /**
         * \brief Constant accessor for the index of the vertex associated with this leaf
         */
        const VertexIndex& vertexIndex() const { return m_vertIndex; }


    public: // Functions related to the associated triangles
        /**
         * \brief Check whether this Leaf has any associated triangles
         */
        bool hasTriangles() const { return !m_tris.empty(); }

        void reserveTriangles(int count) { m_tris.reserve(count); }

        /**
         * \brief Find the number of triangles associated with this leaf
         */
        int numTriangles() const { return m_tris.size(); }



        void addTriangle(TriangleIndex tInd) { m_tris.push_back(tInd); }

        const TriangleList& triangles() const { return m_tris; }
        TriangleList& triangles() { return m_tris;}

    private:
        VertexIndex m_vertIndex;
        TriangleList m_tris;
        bool m_isLeaf;
    };

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
 * \brief Handles generation of the spatial index on a surface mesh
 * for containment queries -- given an arbitrary point in space,
 * is it inside or outside of the surface
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

    typedef Triangle<double, DIM> SpaceTriangle;


    typedef meshtk::Mesh SurfaceMesh;
    typedef meshtk::UnstructuredMesh<meshtk::MIXED> DebugMesh;

    typedef int VertexIndex;
    typedef int TriangleIndex;

    static const int NUM_TRI_VERTS = 3;

    typedef asctoolkit::slam::PositionSet MeshVertexSet;
    typedef asctoolkit::slam::Map<BlockIndex> VertexBlockMap;

    typedef asctoolkit::slam::PositionSet MeshElementSet;
    typedef asctoolkit::slam::Map<SpacePt> VertexPositionMap;

    typedef asctoolkit::slam::policies::CompileTimeStrideHolder<VertexIndex, NUM_TRI_VERTS>  TVStride;
    typedef asctoolkit::slam::StaticConstantRelation<TVStride, MeshElementSet, MeshVertexSet> TriangleVertexRelation;

    typedef typename TriangleVertexRelation::RelationSet TriVertIndices;


    // Some type defs for the Relations from Gray leaf blocks to mesh vertices and elements
    static const int MAX_VERTS_PER_BLOCK = 1;
    typedef asctoolkit::slam::PositionSet GrayLeafSet;
    typedef asctoolkit::slam::policies::CompileTimeStrideHolder<VertexIndex, MAX_VERTS_PER_BLOCK>  BVStride;
    typedef asctoolkit::slam::StaticConstantRelation<BVStride, GrayLeafSet, MeshVertexSet> GrayLeafVertexRelation;
    typedef asctoolkit::slam::StaticVariableRelation                                       GrayLeafElementRelation;
    typedef typename GrayLeafElementRelation::RelationSet TriangleIndexSet;

    typedef asctoolkit::slam::Map<GrayLeafSet>            GrayLeafsLevelMap;
    typedef asctoolkit::slam::Map<GrayLeafVertexRelation> GrayLeafVertexRelationLevelMap;
    typedef asctoolkit::slam::Map<GrayLeafElementRelation> GrayLeafElementRelationLevelMap;


private:
    enum GenerationState { INOUTOCTREE_UNINITIALIZED
                         , INOUTOCTREE_VERTICES_INSERTED
                         , INOUTOCTREE_MESH_REORDERED
                         , INOUTOCTREE_ELEMENTS_INSERTED
                         , INOUTOCTREE_LEAVES_COLORED
    };

public:
    /**
     * \brief Construct an InOutOctree to handle containment queries on a surface mesh
     * \param [in] bb The spatial extent covered by the octree
     */
    InOutOctree(const GeometricBoundingBox& bb, SurfaceMesh*& meshPtr)
        : SpatialOctreeType(bb)
        , m_surfaceMesh(meshPtr)
        , m_vertexSet(0)
        , m_vertexToBlockMap(&m_vertexSet)
        , m_elementSet(0)
        , m_vertexPositions(&m_vertexSet)
        , m_triangleToVertexRelation()
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
     * \param pt The point that we want to check for containment within the surface
     * \return True if the point is within (or on) the surface, false otherwise
     * \note Points outside the octree's bounding box are considered outside the surface
     */
    bool within(const SpacePt& pt) const;

private:

    /**
     * \brief Helper function to insert a vertex into the octree
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


    /**
     * Utility function to print some statistics about the InOutOctree instance
     */
    void printOctreeStats() const;


private:  /** { Mesh Helper functions } */
    /**
     * \brief Helper function to retrieve the position of the vertex from the mesh
     * \param idx The index of the vertex within the surface mesh
     */
    SpacePt getMeshVertexPosition(VertexIndex idx) const
    {
        if(m_surfaceMesh != ATK_NULLPTR)
        {
            SpacePt pt;
            m_surfaceMesh->getMeshNode(idx, pt.data() );
            return pt;
        }
        else
            return m_vertexPositions[idx];
    }

    const SpacePt& vertexPosition(VertexIndex idx) const
    {
        return m_vertexPositions[idx];
    }


    /**
     * Utility function to get the indices of the boundary vertices of a triangle
     * \param idx The index of a triangle within the surface mesh
     */
    TriVertIndices triangleVertexIndices(TriangleIndex idx) const
    {
        return m_triangleToVertexRelation[idx];
    }

    /**
     * \brief Helper function to compute the bounding box of a triangle
     * \param idx The triangle's index within the surface mesh
     */
    GeometricBoundingBox triangleBoundingBox(TriangleIndex idx) const
    {

        GeometricBoundingBox bb;

        // Get the ids of the verts bounding this triangle
        TriVertIndices vertIds = triangleVertexIndices(idx);
        for(int i=0; i< vertIds.size(); ++i)
        {
            bb.addPoint(vertexPosition( vertIds[i] ));
        }

        return bb;
    }

    /**
     * \brief Utility function to retrieve the positions of the traingle's vertices
     * \return A triangle instance whose vertices are positioned in space
     */
    SpaceTriangle trianglePositions(TriangleIndex idx) const
    {
        TriVertIndices verts = triangleVertexIndices(idx);
        return SpaceTriangle( vertexPosition(verts[0])
                            , vertexPosition(verts[1])
                            , vertexPosition(verts[2]));
    }


    VertexIndex leafVertex(const BlockIndex& leafBlk, const InOutBlockData&  leafData) const;

    TriangleIndexSet leafTriangles(const BlockIndex& leafBlk, const InOutBlockData&  leafData) const;


private:

    void dumpOctreeMeshVTK( const std::string& name
                          , bool hasColors = false) const
    {
      #ifdef DUMP_VTK_MESH
        typedef typename OctreeBaseType::MapType LeavesLevelMap;
        typedef typename LeavesLevelMap::const_iterator LeavesIterator;

        DebugMesh* debugMesh= new DebugMesh(3);
        std::stringstream fNameStr;

        int totalLeaves = 0;

        // Iterate through blocks -- count the numbers of leaves
        for(int lev=0; lev< this->m_levels.size(); ++lev)
        {
            const LeavesLevelMap& levelLeafMap = this->m_leavesLevelMap[ lev ];
            for(LeavesIterator it=levelLeafMap.begin(), itEnd = levelLeafMap.end(); it != itEnd; ++it)
            {
                if( it->second.isLeaf())
                    totalLeaves++;
            }
        }

        SLIC_INFO("Dump vtk:: Octree has " << totalLeaves << " leaves.");

        asctoolkit::slam::PositionSet leafSet(totalLeaves);

        typedef asctoolkit::slam::Map<VertexIndex> LeafVertMap;
        typedef asctoolkit::slam::Map<int> LeafTriCountMap;

        LeafVertMap leafVertID(&leafSet);
        LeafVertMap leafVertID_unique(&leafSet);
        LeafTriCountMap leafTriCount(&leafSet);
        LeafTriCountMap leafColors(&leafSet);

        // Iterate through blocks -- count the numbers of leaves
        int leafCount = 0;
        for(int lev=0; lev< this->m_levels.size(); ++lev)
        {
            const LeavesLevelMap& levelLeafMap = this->m_leavesLevelMap[ lev ];
            for(LeavesIterator it=levelLeafMap.begin(), itEnd = levelLeafMap.end(); it != itEnd; ++it)
            {
                if( ! it->second.isLeaf())
                    continue;


                // Add a hex
                BlockIndex block( it->first, lev);
                const InOutLeafData& leafData = it->second;
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

                debugMesh->insertCell( data, meshtk::LINEAR_HEX, 8);

                int vIdx = leafData.vertexIndex();
                SpacePt pt = vertexPosition(vIdx);

                //SLIC_INFO("In block " << block<< " vertex id is " << vIdx
                //      << " associated point is " << pt);


                leafVertID[leafCount] = vIdx;
                leafVertID_unique[leafCount] = (leafData.hasVertex() && blockIndexesVertex(vIdx, block))
                        ? vIdx
                        : InOutLeafData::NO_VERTEX;
                leafTriCount[leafCount] = leafData.numTriangles();

                if(hasColors)
                {
                    switch( leafData.color() )
                    {
                    case InOutLeafData::Black:          leafColors[leafCount] =  1; break;
                    case InOutLeafData::White:          leafColors[leafCount] = -1; break;
                    case InOutLeafData::Gray:           leafColors[leafCount] =  0; break;
                    case InOutLeafData::Undetermined:
                        SLIC_ASSERT_MSG(false, "Leaf " << block << " of InOutOctree is missing a color.");
                        break;
                    }
                }

                leafCount++;
            }
        }

        // Add the fields
        meshtk::FieldData* CD = debugMesh->getCellFieldData();
        CD->addField( new meshtk::FieldVariable< VertexIndex >("vertID", leafSet.size()) );
        CD->addField( new meshtk::FieldVariable< VertexIndex >("uniqVertID", leafSet.size()) );
        CD->addField( new meshtk::FieldVariable< int >("triCount", leafSet.size()) );

        if(hasColors)
            CD->addField( new meshtk::FieldVariable< int >("colors", leafSet.size()) );

        VertexIndex* vertID = CD->getField( "vertID" )->getIntPtr();
        VertexIndex* uniqVertID = CD->getField( "uniqVertID" )->getIntPtr();
        int* triCount = CD->getField( "triCount" )->getIntPtr();
        int* colors = hasColors ? CD->getField( "colors" )->getIntPtr() : ATK_NULLPTR;

        SLIC_ASSERT( vertID != ATK_NULLPTR );
        SLIC_ASSERT( uniqVertID != ATK_NULLPTR );
        SLIC_ASSERT( triCount != ATK_NULLPTR );
        SLIC_ASSERT( !hasColors || (colors != ATK_NULLPTR) );

        for ( int i=0; i < leafSet.size(); ++i ) {
            vertID[i] = leafVertID[i];
            uniqVertID[i] = leafVertID_unique[i];
            triCount[i] = leafTriCount[i];

            if(hasColors)
                colors[i] = leafColors[i];
        }

        debugMesh->toVtkFile(name);

        delete debugMesh;
        debugMesh = ATK_NULLPTR;
      #else
        // Do something with the parameters to avoid a warning about unused parameters
        ATK_DEBUG_VAR(name);
        ATK_DEBUG_VAR(hasColors);
      #endif
    }

    /**
     * \brief Utility function to dump a single element vtk file
     */
    void dumpMeshVTK( const std::string& name
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
            SpaceTriangle triPos = trianglePositions(idx);

            debugMesh->insertNode( triPos.A()[0], triPos.A()[1], triPos.A()[2]);
            debugMesh->insertNode( triPos.B()[0], triPos.B()[1], triPos.B()[2]);
            debugMesh->insertNode( triPos.C()[0], triPos.C()[1], triPos.C()[2]);

            int verts[3] = {0,1,2};
            debugMesh->insertCell(verts, meshtk::LINEAR_TRIANGLE, 3);

            fNameStr << name << idx << ".vtk";

            SLIC_INFO("// Triangle " << idx );
            SLIC_INFO("TriangleType tri(PointType::make_point" << triPos.A()
                                  << ", PointType::make_point" << triPos.B()
                                  << ", PointType::make_point" << triPos.C()<<");" );

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

            debugMesh->insertCell( data, meshtk::LINEAR_HEX, 8);

            fNameStr << name << idx << "_" << counter[name] << ".vtk";

            SLIC_INFO("// Block index " << block);
            SLIC_INFO("BoundingBoxType box"<<counter[name]<<"(PointType::make_point" << blockBB.getMin() << ", PointType::make_point" << blockBB.getMax()<<");" );
            counter[name]++;
        }

        debugMesh->toVtkFile(fNameStr.str());

        delete debugMesh;
        debugMesh = ATK_NULLPTR;
      #else
        // Do something with the parameters to avoid a warning about unused parameters
        ATK_DEBUG_VAR(name);
        ATK_DEBUG_VAR(idx);
        ATK_DEBUG_VAR(block);
        ATK_DEBUG_VAR(blockBB);
        ATK_DEBUG_VAR(isTri);
      #endif
    }

    /**
     * Deal with numerical degeneracies on triangle boundaries by explicitly checking for
     * block containing triangle vertices
     */
    bool blockIndexesVertex(TriangleIndex idx, const BlockIndex& block) const
    {
        TriVertIndices tVerts = triangleVertexIndices(idx);

        if(m_generationState >= INOUTOCTREE_MESH_REORDERED)
        {
#if 0
            const InOutBlockData& leafData = (*this)[block];
            return ( leafData.hasData()) ? incidentInVertex(tVerts, leafVertex(block,leafData)) : false;
#else
            for(int i=0; i< tVerts.size(); ++i)
            {
                if( block == m_vertexToBlockMap[ tVerts[i]])
                    return true;
            }
#endif
        }
        else
        {
            for(int i=0; i< tVerts.size(); ++i)
            {
                if( block == this->findLeafBlock( vertexPosition(tVerts[i])))
                    return true;
            }
        }

        return false;
    }

    /**
     * \brief Checks whether the indexed triangle contains a reference to the given vertex
     */
    bool incidentInVertex(const TriVertIndices& triVerts, VertexIndex vIdx) const
    {
        return (triVerts[0] == vIdx) || (triVerts[1] == vIdx) || (triVerts[2] == vIdx);
    }


private:
    /**
     * \brief Checks if all indexed triangles in the block share a common vertex
     * \param leafBlock [in] The current octree block
     * \param leafData [inout] The data associated with this block
     * \note A side effect of this function is that we set the leafData's vertex to the common
     * vertex if one is found
     * \return True, if all triangles indexed by this leaf share a common vertex, false otherwise.
     */
    bool allTrianglesIncidentInCommonVertex(const BlockIndex& leafBlock, DynamicGrayBlockData& leafData) const;

    /**
     * \brief Finds a color for the given block blk and propagates to neighbors
     * \note Propagates color to same-level and coarser level neighbors
     * \param blk The block to color
     * \param blkData The data associated with this block
     * \return True if we were able to find a color for blk, false otherwise
     */
    bool colorLeafAndNeighbors(const BlockIndex& blk, InOutBlockData& blkData);


private:
    DISABLE_COPY_AND_ASSIGNMENT(InOutOctree);

    /** \brief Checks internal consistency of the octree representation */
    void checkValid() const;


    /**
     * \brief Determines whether the specified point is within the gray leaf
     * \param pt The point we are querying
     * \param leafBlk The block of the gray leaf
     * \param data The data associated with the leaf block
     * \return True, if the point is inside the local surface associated with this block, false otherwise
     */
    bool withinGrayBlock( const SpacePt & pt, const BlockIndex& leafBlk, const InOutBlockData&  data) const;

  /** \brief Add the vertex positions and triangle boundary relations to the surface mesh */
    void regenerateSurfaceMesh();


protected:
    SurfaceMesh*& m_surfaceMesh; /** pointer to pointer to allow changing the mesh */

    MeshVertexSet m_vertexSet;
    VertexBlockMap m_vertexToBlockMap;

    MeshElementSet m_elementSet;
    VertexPositionMap m_vertexPositions;
    TriangleVertexRelation m_triangleToVertexRelation;

    GrayLeafsLevelMap m_grayLeafsMap;
    GrayLeafVertexRelationLevelMap m_grayLeafToVertexRelationLevelMap;
    GrayLeafElementRelationLevelMap m_grayLeafToElementRelationLevelMap;

    GenerationState m_generationState;
};


namespace{
#ifdef ATK_DEBUG
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
    typedef asctoolkit::utilities::Timer Timer;

    // Loop through mesh vertices
    SLIC_INFO("Generating index on mesh vertices.");

    Timer timer;

    // STEP 1 -- Add mesh vertices to octree
    timer.start();
    int numMeshVerts = m_surfaceMesh->getMeshNumberOfNodes();
    for(int idx=0; idx < numMeshVerts; ++idx)
    {
        insertVertex(idx);
    }
    timer.stop();
    m_generationState = INOUTOCTREE_VERTICES_INSERTED;
    SLIC_INFO("\tInserting vertices took " << timer.elapsed() << " seconds.");

    // STEP 1(b) -- Update the mesh vertices and cells with after vertex welding from octree
    timer.start();
    updateSurfaceMeshVertices();
    timer.stop();
    m_generationState = INOUTOCTREE_MESH_REORDERED;
    SLIC_INFO("\tUpdating mesh took " << timer.elapsed() << " seconds.");


  #ifdef DUMP_OCTREE_INFO
    // -- Print some stats about the octree
    SLIC_INFO("** Octree stats after inserting vertices");
    dumpOctreeMeshVTK("prOctree.vtk", false);
    printOctreeStats();
  #endif
    checkValid();


    // STEP 2 -- Add mesh triangles to octree
    timer.start();
    insertMeshTriangles();
    timer.stop();
    m_generationState = INOUTOCTREE_ELEMENTS_INSERTED;
    SLIC_INFO("\tInserting triangles took " << timer.elapsed() << " seconds.");

    // STEP 3 -- Color the blocks of the octree -- Black (in), White(out), Gray(Intersects surface)
    timer.start();
    colorOctreeLeaves();

    timer.stop();
    m_generationState = INOUTOCTREE_LEAVES_COLORED;
    SLIC_INFO("\tColoring octree leaves took " << timer.elapsed() << " seconds.");

    // -- Print some stats about the octree
  #ifdef DUMP_OCTREE_INFO
    SLIC_INFO("** Octree stats after inserting triangles");
    dumpOctreeMeshVTK("pmOctree.vtk", true);
    printOctreeStats();
  #endif
    checkValid();

    // CLEANUP -- Finally, fix up the surface mesh after octree operations
    SLIC_INFO("\tRegenerating the mesh");
    regenerateSurfaceMesh();

    SLIC_INFO("\tFinished generating the InOutOctree.");
}


template<int DIM>
void InOutOctree<DIM>::insertVertex (VertexIndex idx, int startingLevel)
{
    static const double EPS_SQ = 1e-18;


    const SpacePt pt = getMeshVertexPosition(idx);

    BlockIndex block = this->findLeafBlock(pt, startingLevel);
    InOutBlockData& blkData = (*this)[block];


    #ifdef ATK_DEBUG
    if(idx == DEBUG_VERT_IDX)
    {
        SLIC_DEBUG("\t -- inserting pt with index " << idx << " and coords: " << pt
                  << ". Looking at block " << block
                  << " w/ blockBB " << this->blockBoundingBox(block)
                  << " indexing leaf vertex " << blkData.dataIndex()
                    );
    }
    #endif


    if(! blkData.hasData() )
    {
        blkData.setData(idx);

        // Update the vertex-to-block map for this vertex
        if(!m_vertexSet.empty())
            m_vertexToBlockMap[idx] = block;
    }
    else
    {
        // check if we should merge the vertices
        VertexIndex origVertInd = blkData.dataIndex();
        if( squared_distance( pt, getMeshVertexPosition(origVertInd) ) >= EPS_SQ )
        {
            blkData.clear();
            this->refineLeaf(block);

            insertVertex(origVertInd, block.childLevel() );
            insertVertex(idx, block.childLevel() );
        }
    }

    #ifdef ATK_DEBUG
    if( blkData.dataIndex() == DEBUG_VERT_IDX)
    {
        SLIC_DEBUG("-- vertex " << idx
                  << " is indexed in block " << block
                  << ". Leaf vertex is " << blkData.dataIndex() );
    }
    #endif
}


template<int DIM>
void InOutOctree<DIM>::insertMeshTriangles ()
{
    typedef asctoolkit::utilities::Timer Timer;

    typedef typename OctreeBaseType::MapType LeavesLevelMap;
    typedef typename LeavesLevelMap::iterator LeavesIterator;

    SLIC_ASSERT( !m_vertexSet.empty() );

    // Temporary usage of DyamicGrayBlockData for current level
//#if defined(USE_CXX11)
//    typedef std::unordered_map<int, DynamicGrayBlockData> DynamicLevelData;
//#else
//    typedef boost::unordered_map<int, DynamicGrayBlockData> DynamicLevelData;
//#endif
    typedef std::vector<DynamicGrayBlockData> DynamicLevelData;



    DynamicLevelData currentLevelData;
    DynamicLevelData nextLevelData;
    currentLevelData.reserve( 1<<10);
    nextLevelData.reserve( 1<<10);

    std::vector<VertexIndex>   gvRelData;
    std::vector<int>           geSizeRelData;
    std::vector<TriangleIndex> geIndRelData;


    /// --- Initialize root level data
    // Add all triangle references to the root
    BlockIndex rootBlock = this->root();
    InOutBlockData& rootBlkData = (*this)[rootBlock];

    DynamicGrayBlockData rootData;
    rootData.triangles().reserve( m_elementSet.size() );
    for(int idx=0; idx < m_elementSet.size(); ++idx)
        rootData.addTriangle(idx);

    if( rootBlkData.hasData() )
        rootData.setVertex( rootBlkData.dataIndex() );
    rootData.setLeafFlag(rootBlkData.isLeaf());

    // Update tree and current level map with new information
    rootBlkData.setData(0);
    currentLevelData.push_back(rootData);


    /// --- Iterate through octree levels and insert triangles into the blocks that they intersect
    for(int lev=0; lev< this->m_levels.size(); ++lev)
    {
        Timer levelTimer(true);

        gvRelData.clear();
        geSizeRelData.clear();
        geSizeRelData.push_back(0);

        geIndRelData.clear();

        int nextLevelDataBlockCounter = 0;

        LeavesLevelMap& levelLeafMap = this->m_leavesLevelMap[ lev ];

        for(LeavesIterator it=levelLeafMap.begin(), itEnd = levelLeafMap.end(); it != itEnd; ++it)
        {
            InOutBlockData& blkData = it->second;

            if(! blkData.hasData() )
                continue;

            BlockIndex blk(it->first, lev);
            DynamicGrayBlockData& dynamicLeafData = currentLevelData[ blkData.dataIndex() ];


            // Get the associated vertex -- for leaf blocks; does not exist for internal blocks

            // If leaf
            //  -- check if we must refine -- if so, refine and add vertex to child
            //  -- else, increment gray counter and add the associated data to the relations,


            //SLIC_INFO("Adding triangles to blk " << blk << " with data " << blkData.dataIndex() );

            bool isInternal = !dynamicLeafData.isLeaf();

            //if(isInternal)
            //    SLIC_INFO("\tblock is internal");


            bool isLeafThatMustRefine = !isInternal && ! allTrianglesIncidentInCommonVertex(blk, dynamicLeafData);

            //if(isLeafThatMustRefine)
            //    SLIC_INFO("\bblock must refine");


            // Leaf blocks that don't refine are 'finalized' -- add  them to the relation
            if( !isInternal && !isLeafThatMustRefine )
            {
                if( dynamicLeafData.hasTriangles())
                {
//                    SLIC_INFO("\tBlock " << blk << " is a leaf block that is not refining -- "
//                            << " adding as " << gvRelData.size() << "th gray node at level " << lev
//                            << ". blk has " << dynamicLeafData.numTriangles() << " triangles."
//                            );

                    blkData.setData( gvRelData.size() );

                    gvRelData.push_back( dynamicLeafData.vertexIndex() );

                    std::copy( dynamicLeafData.triangles().begin()
                             , dynamicLeafData.triangles().end()
                             , std::back_inserter(geIndRelData));

                    geSizeRelData.push_back( geIndRelData.size());
                }
            }
            else
            {
                // At this point, we must distribute the triangles among the appropriate children
                DynamicGrayBlockData::TriangleList& parentTriangles = dynamicLeafData.triangles();
                int numTriangles = parentTriangles.size();

                // Refine the leaf if necessary
                if( isLeafThatMustRefine )
                {
                    VertexIndex vIdx = dynamicLeafData.vertexIndex();

//                    SLIC_INFO("Refining leaf " << blk << " with vertex " << vIdx);

                    this->refineLeaf(blk);
                    dynamicLeafData.setLeafFlag(false);

                    // Note: At this stage, leaf might have triangles but no assigned vertex,
                    //       so we need to check for that
                    if( vIdx >= 0 && m_vertexToBlockMap[vIdx] == blk )
                        insertVertex(vIdx, blk.childLevel());
                }
                else if(isInternal)
                {
                    blkData.setInternal();
                }

                // Cache information about children blocks
                BlockIndex              childBlk[BlockIndex::NUM_CHILDREN];
                GeometricBoundingBox    childBB[BlockIndex::NUM_CHILDREN];
                DynamicGrayBlockData    childData[BlockIndex::NUM_CHILDREN];
                DynamicGrayBlockData*   childDataPtr[BlockIndex::NUM_CHILDREN];

                for(int j=0; j< BlockIndex::NUM_CHILDREN; ++j)
                {
                    childBlk[j]  = blk.child(j);

                    SLIC_ASSERT_MSG( this->hasBlock( childBlk[j] )
                                   , "Octree does not have child " << j << " -- " << childBlk[j]
                                   << "-- of block " << blk);



                    childBB[j]   = this->blockBoundingBox( childBlk[j] );

                    InOutBlockData& iob = (*this)[childBlk[j] ];

                    if(!iob.hasData())
                    {
                        childData[j] = DynamicGrayBlockData();
                        childData[j].setLeafFlag( iob.isLeaf() );
                    }
                    else
                    {
                        childData[j] = DynamicGrayBlockData(iob.dataIndex(), iob.isLeaf());
                    }

                    childDataPtr[j] = & childData[j];
                }

                if(nextLevelData.capacity() < ( nextLevelData.size() + BlockIndex::NUM_CHILDREN ))
                    nextLevelData.reserve(nextLevelData.size() * 4);

                // Add all triangles to intersecting children blocks
                for(int i=0; i< numTriangles; ++i)
                {
                    TriangleIndex tIdx = parentTriangles[i];
                    SpaceTriangle spaceTri = trianglePositions(tIdx);
                    GeometricBoundingBox tBB = triangleBoundingBox(tIdx);
                    //tBB.scale(1.00001);

                    for(int j=0; j< BlockIndex::numChildren(); ++j)
                    {
                        bool shouldAddTriangle = childDataPtr[j]->isLeaf()
                                ? blockIndexesVertex(tIdx, childBlk[j]) || intersect(spaceTri, childBB[j])
                                : intersect(tBB, childBB[j])
                                ;

                        if(shouldAddTriangle)
                        {
                            // Place the DynamicGrayBlockData in the array before adding its data
                            if(!childDataPtr[j]->hasTriangles())
                            {
                                // Copy the DynamicGrayBlockData into the array
                                nextLevelData.push_back(childData[j]);

                                // Update the pointer
                                childDataPtr[j] = &nextLevelData[nextLevelDataBlockCounter];

                                // Set the data in the octree to this index and update the index
                                (*this)[childBlk[j]].setData(nextLevelDataBlockCounter++);
                            }

                            childDataPtr[j]->addTriangle(tIdx);
                        }
                    }
                }
            }
        }

        if(! levelLeafMap.empty() )
        {
//            SLIC_DEBUG("\t -- About to create gray leaf set of size " << gvRelData.size() << " for level " << lev);
//            SLIC_DEBUG("\t -- size of grayleaf to vertex relation " << gvRelData.size());
//            SLIC_DEBUG("\t -- size of grayleaf to element sizes " << geSizeRelData.size());
//            SLIC_DEBUG("\t -- size of grayleaf to element ind array " << geIndRelData.size());

            m_grayLeafsMap[lev] = GrayLeafSet( gvRelData.size() );
            m_grayLeafToVertexRelationLevelMap[lev] = GrayLeafVertexRelation(&m_grayLeafsMap[lev], &m_vertexSet);
            m_grayLeafToVertexRelationLevelMap[lev].bindRelationData(gvRelData);

            m_grayLeafToElementRelationLevelMap[lev] = GrayLeafElementRelation(&m_grayLeafsMap[lev], &m_elementSet);
            m_grayLeafToElementRelationLevelMap[lev].bindRelationData(geSizeRelData, geIndRelData);
        }

        currentLevelData.clear();
        nextLevelData.swap(currentLevelData);

        if(! levelLeafMap.empty() )
          SLIC_DEBUG("\tInserting triangles into level " << lev <<" took " << levelTimer.elapsed() <<" seconds.");
    }
}

template<int DIM>
void InOutOctree<DIM>::regenerateSurfaceMesh()
{
    if(m_surfaceMesh != ATK_NULLPTR)
    {
        delete m_surfaceMesh;
        m_surfaceMesh = ATK_NULLPTR;
    }

    typedef meshtk::UnstructuredMesh< meshtk::LINEAR_TRIANGLE > TriangleMesh;
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
        triMesh->insertCell(tv, meshtk::LINEAR_TRIANGLE, 3);
    }

    m_surfaceMesh = triMesh;
}

template<int DIM>
bool InOutOctree<DIM>::colorLeafAndNeighbors(const BlockIndex& leafBlk, InOutBlockData& leafData)
{
    bool isColored = leafData.isColored();

    if( ! isColored )
    {
        // Leaf does not yet have a color... try to find its color from same-level face neighbors
        for(int i=0; !isColored && i< leafBlk.numFaceNeighbors(); ++i)
        {
            BlockIndex neighborBlk = leafBlk.faceNeighbor(i);
            if( this->isLeaf( neighborBlk ) )
            {
                const InOutBlockData& neighborData = (*this)[neighborBlk];
                switch(neighborData.color())
                {
                case InOutBlockData::Black: leafData.setBlack(); break;
                case InOutBlockData::White: leafData.setWhite(); break;
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
                    switch(leafData.color())
                    {
                    case InOutBlockData::Black:
                        neighborData.setBlack();
                        break;
                    case InOutBlockData::White:
                        neighborData.setWhite();
                        break;
                    case InOutBlockData::Gray:
                        // PROBLEM HERE -- some neighbors are getting the wrong color!
                        if( withinGrayBlock(this->blockBoundingBox(neighborBlk).centroid(),  leafBlk, leafData) )
                            neighborData.setBlack();
                        else
                            neighborData.setWhite();
                        break;
                    case InOutBlockData::Undetermined:
                        break;
                    }
                }
            }
        }
    }

    return isColored;
}

template<int DIM>
void InOutOctree<DIM>::colorOctreeLeaves()
{
    typedef asctoolkit::utilities::Timer Timer;

    typedef typename OctreeBaseType::MapType LeavesLevelMap;
    typedef typename LeavesLevelMap::iterator LeavesIterator;

    // Bottom-up traversal of octree
    for(int lev=this->maxLeafLevel()-1; lev >= 0; --lev)
    {
        Timer levelTimer(true);

        typedef std::vector<GridPt> GridPtVec;
        GridPtVec uncoloredBlocks;

        LeavesLevelMap& levelLeafMap = this->m_leavesLevelMap[ lev ];
        for(LeavesIterator it=levelLeafMap.begin(), itEnd = levelLeafMap.end(); it != itEnd; ++it)
        {
            if( !it->second.isLeaf() )
                continue;


            BlockIndex leafBlk(it->first, lev);

            if(! colorLeafAndNeighbors( leafBlk, it->second) )
                uncoloredBlocks.push_back( it->first);

//            // Verified assumption: Existence of leaf implies that either
//            // * it is gray
//            // * one of its siblings is gray
//            // * one of its siblings has a gray descendant


        }

        // Iterate through the uncolored blocks until all have a color
        // This should terminate since we know that one of its siblings (or their descendants) is gray
        while(! uncoloredBlocks.empty())
        {
            //SLIC_INFO("\tAt level " << lev << ": Attempting to add " << uncoloredBlocks.size() <<  " uncolored block");


            int prevCount = uncoloredBlocks.size();
            ATK_DEBUG_VAR(prevCount);

            GridPtVec prevVec;
            prevVec.swap(uncoloredBlocks);
            for(typename GridPtVec::iterator it = prevVec.begin(), itEnd = prevVec.end(); it < itEnd; ++it)
            {
                BlockIndex leafBlk(*it, lev);
                if(! colorLeafAndNeighbors( leafBlk, (*this)[leafBlk] ) )
                    uncoloredBlocks.push_back( *it);
            }

            SLIC_ASSERT( static_cast<int>(uncoloredBlocks.size()) < prevCount);
        }


      #ifdef ATK_DEBUG
        ////Assumption is that by the time we get to the end, all blocks will be colored...
        for(LeavesIterator it=levelLeafMap.begin(), itEnd = levelLeafMap.end(); it != itEnd; ++it)
        {
            if( ! it->second.isLeaf() )
                continue;

            SLIC_ASSERT_MSG( it->second.isColored()
                            , "Error after coloring level " << lev
                            << " leaf block " << BlockIndex(it->first, lev)
                            << " was not colored."
            );
        }
       #endif

        if(! levelLeafMap.empty() )
            SLIC_DEBUG("\tColoring level "<< lev << " took " << levelTimer.elapsed() << " seconds.");

    }
}


template<int DIM>
typename InOutOctree<DIM>::VertexIndex InOutOctree<DIM>::leafVertex(const BlockIndex& leafBlk, const InOutBlockData&  leafData) const
{
//    SLIC_DEBUG("Inside leafVertex -- with blk " << leafBlk <<" and data " << leafData.dataIndex());

    if( m_generationState >= INOUTOCTREE_ELEMENTS_INSERTED)
    {
        SLIC_ASSERT( m_grayLeafToVertexRelationLevelMap[leafBlk.level()].size() > 0 );
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
//    SLIC_ASSERT( m_grayLeafToElementRelationLevelMap[leafBlk.level()].size() > 0 );

    return m_grayLeafToElementRelationLevelMap[leafBlk.level()][leafData.dataIndex()];
}



template<int DIM>
bool InOutOctree<DIM>::withinGrayBlock(const SpacePt & pt, const BlockIndex& leafBlk, const InOutBlockData&  leafData) const
{
    SLIC_ASSERT( leafData.color() == InOutBlockData::Gray );
    SLIC_ASSERT( leafData.hasData() );

    VertexIndex vIdx = leafVertex(leafBlk, leafData);

    //SLIC_INFO("withinGray -- vertex index is " << vIdx << " for leaf block " << leafBlk << " with data " << leafData.dataIndex());


    const SpaceVector vec(vertexPosition( vIdx), pt );

    // Iterate through triangles, count signs of dot products of vec with normals
    // If all dot products are positive, point is outside
    // If all dot products are negative, point is inside
    // Otherwise, short-circuit the calculation
    int sgnCount = 0;

    TriangleIndexSet triSet = leafTriangles(leafBlk, leafData);
    const int numTris = triSet.size();
    for(int i=0; i< numTris; ++i)
    {
        TriangleIndex tIdx = triSet[i];

        if( vec.dot( trianglePositions(tIdx).normal() ) >= 0 )
        {
            if(sgnCount == i)
                ++sgnCount;
            else
                break;
        }
        else
        {
            if(sgnCount == -i)
                --sgnCount;
            else
                break;
        }
    }

    if(sgnCount == numTris)             // outside
        return false;
    else if(sgnCount == -numTris)       // inside
        return true;


    // Else, pt is within some half-spaces and outside others
    // We determine the containment of the point w.r.t. the surface
    // based on whether the geometry in the gray block is locally convex
    //  and point is inside if surface is locally concave (inside) or convex (outside)
    // We can determine this via an orientation test on two of the triangles.

    SLIC_ASSERT( numTris > 1 );

    TriangleIndex firstTriangle = triSet[0];

    // Find a vertex from the local surface that is not incident in the first triangle
    TriVertIndices firstTriVerts = triangleVertexIndices( firstTriangle );
    TriVertIndices lastTriVerts = triangleVertexIndices( triSet[ numTris -1]);
    VertexIndex otherVert = lastTriVerts[0];
    if( incidentInVertex(firstTriVerts, otherVert))
        otherVert = lastTriVerts[1];
    if( incidentInVertex(firstTriVerts, otherVert))
        otherVert = lastTriVerts[2];
    SLIC_ASSERT( !incidentInVertex(firstTriVerts, otherVert) );

    return ON_NEGATIVE_SIDE != orientation( vertexPosition(otherVert), trianglePositions( firstTriangle));
}







template<int DIM>
void InOutOctree<DIM>::updateSurfaceMeshVertices()
{
    typedef asctoolkit::slam::Map<VertexIndex> IndexMap;

    // Create a map from old vertex indices to new vertex indices
    MeshVertexSet oldVerts( m_surfaceMesh->getMeshNumberOfNodes() );
    const int NO_DATA = -1;
    IndexMap vertexIndexMap( &oldVerts, NO_DATA);

    // Generate unique indices for new mesh vertices
    int uniqueVertexCounter = 0;
    for(int i=0; i< oldVerts.size(); ++i)
    {
        // Find the block and its indexed vertex in the octree
        BlockIndex leafBlock = this->findLeafBlock( getMeshVertexPosition(i) );
        SLIC_ASSERT( (*this)[leafBlock].hasData() );
        VertexIndex vInd = (*this)[ leafBlock ].dataIndex();

        // If the indexed vertex doesn't have a new id, give it one
        if(vertexIndexMap[vInd] == NO_DATA)
        {
            vertexIndexMap[vInd] = uniqueVertexCounter++;
        }

        // If this is not the indexed vertex, grab that vertex's new IDX
        if(vInd != i)
        {
            vertexIndexMap[i] = vertexIndexMap[vInd];
        }
    }

    // Create a vertex set on the new vertices and grab coordinates from the old ones
    m_vertexSet = MeshVertexSet( uniqueVertexCounter );
    m_vertexPositions = VertexPositionMap(&m_vertexSet);
    for(int i=0; i< oldVerts.size(); ++i)
    {
        const VertexIndex& vInd = vertexIndexMap[i];
        m_vertexPositions[vInd] = getMeshVertexPosition(i);
    }

    // Update the octree leaf vertex IDs to the new mesh
    // and create the map from the new vertices to their octree blocks
    m_vertexToBlockMap = VertexBlockMap(&m_vertexSet);
    for(int i = 0; i< m_vertexSet.size(); ++i)
    {
        const SpacePt& pos = m_vertexPositions[i];
        BlockIndex leafBlock = this->findLeafBlock(pos);
        SLIC_ASSERT( this->isLeaf(leafBlock) && (*this)[leafBlock].hasData() );

        (*this)[ leafBlock ].setData(i);
        m_vertexToBlockMap[i] = leafBlock;
    }

    // Update the vertex IDs of the triangles to the new vertices and create a SLAM relation on these
    int numOrigTris = m_surfaceMesh->getMeshNumberOfCells();
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
    m_surfaceMesh = ATK_NULLPTR;


    SLIC_INFO("After inserting vertices, mesh has "
              << m_vertexSet.size() << " vertices" << " and "
              << m_elementSet.size() << " triangles." );

}

template<int DIM>
bool InOutOctree<DIM>::allTrianglesIncidentInCommonVertex(const BlockIndex& leafBlock
                                                         , DynamicGrayBlockData& leafData) const
{
    bool shareCommonVert = false;
    static const VertexIndex NO_VERTEX = -1;

    VertexIndex commonVert = NO_VERTEX;

    const int numTris = leafData.numTriangles();
    const DynamicGrayBlockData::TriangleList& tris = leafData.triangles();

    if(numTris == 0)
    {
        shareCommonVert = true;     // No need to refine here
    }
    else if(leafData.hasVertex() && m_vertexToBlockMap[ leafData.vertexIndex()] == leafBlock )
    {
        // This is a leaf node containing a vertex
        // Loop through the triangles and check that all are incident with this vertex
        commonVert = leafData.vertexIndex();
        shareCommonVert = true;
        for(int i=0; shareCommonVert && i< numTris; ++i)
        {
            if( !incidentInVertex(triangleVertexIndices(tris[i]), commonVert) )
            {
                shareCommonVert = false;
            }
        }
    }
    else
    {
        // Simple strategy -- can make more efficient later
        // Note: Assumption is that we have very few (i.e. 6 or fewer) triangles in a bucket
        //       Optimize for cases where there are only 1 or 2 indexed triangles
        if(numTris == 1)
        {
            commonVert = triangleVertexIndices(tris[0])[0];
            shareCommonVert = true;
        }
        else if(numTris == 2)
        {
            // There are two triangles -- check that they have at least one common vertex
            TriVertIndices tvRel0 = triangleVertexIndices(tris[0]);
            TriVertIndices tvRel1 = triangleVertexIndices(tris[1]);

            for(int i=0; !shareCommonVert && i< tvRel1.size(); ++i)
            {
                commonVert = tvRel1[i];
                if( incidentInVertex(tvRel0, commonVert) )
                    shareCommonVert = true;
            }
        }
        else    // numTris >= 3
        {
            /// Step 1: Find a vertex that the first three triangles share
            int sharedVert = NO_VERTEX;
            std::vector<VertexIndex> t01Verts;
            t01Verts.reserve(2);

            // Compare the first two triangles -- there can be at most two matches
            TriVertIndices t0Verts = triangleVertexIndices(tris[0]);
            TriVertIndices t1Verts = triangleVertexIndices(tris[1]);
            for(int i=0; i< t1Verts.size(); ++i)
            {
                if( incidentInVertex( t0Verts, t1Verts[i]))
                    t01Verts.push_back(t1Verts[i]);
            }

            // Now check the third triangle's vertices against these matches
            TriVertIndices t2Verts = triangleVertexIndices(tris[2]);
            for(std::size_t i=0; !shareCommonVert && i< t01Verts.size(); ++i)
            {
                if( incidentInVertex( t2Verts, t01Verts[i]))
                {
                    sharedVert = t01Verts[i];
                    shareCommonVert = true;
                }
            }

            /// Step 2: If we got to here, check that all other triangles have this vertex
            for(int i=3; shareCommonVert && i< numTris; ++i)
            {
                if( ! incidentInVertex( triangleVertexIndices(tris[i]), sharedVert) )
                        shareCommonVert = false;
            }

            if(shareCommonVert)
            {
                commonVert = sharedVert;
            }
        }
    }

    leafData.setVertex( commonVert );
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
        case InOutBlockData::Black:      return true;
        case InOutBlockData::White:      return false;
        case InOutBlockData::Gray:       return withinGrayBlock( pt, block, data);
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
    typedef typename OctreeBaseType::MapType LeavesLevelMap;
    typedef typename LeavesLevelMap::const_iterator LeavesIterator;


    typedef asctoolkit::slam::Map<int> LeafCountMap;

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
    //
    for(int lev=0; lev< this->m_levels.size(); ++lev)
    {
        const LeavesLevelMap& levelLeafMap = this->m_leavesLevelMap[ lev ];
        levelBlocks[lev] = levelLeafMap.size();
        levelLeaves[lev] = 0;
        levelLeavesWithVert[lev] = 0;
        levelTriangleRefCount[lev] = 0;
        levelWhiteBlockCount[lev] = 0;
        levelBlackBlockCount[lev] = 0;
        levelGrayBlockCount[lev] = 0;

        for(LeavesIterator it=levelLeafMap.begin(), itEnd = levelLeafMap.end(); it != itEnd; ++it)
        {
            const InOutBlockData& blockData = it->second;
            BlockIndex block(it->first, lev);

            if(blockData.isLeaf())
            {
                ++levelLeaves[lev];

                if(blockData.hasData())
                {
                    ++levelLeavesWithVert[ lev ];

                    if(m_generationState >= INOUTOCTREE_ELEMENTS_INSERTED)
                    {
                        TriangleIndexSet triSet = leafTriangles(block, blockData);
                        if(triSet.size() >= 0)
                            levelTriangleRefCount[ lev ] += triSet.size();
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

            octreeStatsStr <<"\n";
        }
    }


    double meshNumTriangles = static_cast<double>(m_elementSet.size());
    int percentWithVert = (100. * totalLeavesWithVert) / totalLeaves;
    octreeStatsStr<<"  Mesh has " << m_vertexSet.size() << " vertices."
                 <<"\n  Octree has " << totalBlocks << " blocks; "
                 <<  totalBlocks - totalLeaves << " internal; "
                 <<  totalLeaves << " leaves "
                 << " (" <<   percentWithVert  << "% w/ vert); ";
    if(m_generationState >= INOUTOCTREE_ELEMENTS_INSERTED)
    {
        octreeStatsStr<<" \n\t There were " << totalTriangleRefCount << " triangle references "
                 <<" (avg. " << ( totalTriangleRefCount /  (double)meshNumTriangles ) << " refs per triangle).";
    }
    if(m_generationState >= INOUTOCTREE_LEAVES_COLORED)
    {
        octreeStatsStr<<"\n Colors B,W,G ==> " << totalBlackBlocks
                 << "," << totalWhiteBlocks
                 << "," << totalGrayBlocks
                 <<"\n";
    }

    SLIC_INFO( octreeStatsStr.str() );

    if(m_generationState >= INOUTOCTREE_LEAVES_COLORED)
    {
        if(DEBUG_TRI_IDX >= 0)
        {
             dumpMeshVTK("triangle", DEBUG_TRI_IDX, this->root(), this->blockBoundingBox(this->root()), true);


             TriVertIndices tv = triangleVertexIndices(DEBUG_TRI_IDX);
             for(int i=0; i< 3; ++i)
             {
                 BlockIndex blk = m_vertexToBlockMap[ tv[i] ];
                 dumpMeshVTK("triangleVertexBlock", DEBUG_TRI_IDX, blk, this->blockBoundingBox(blk), false);
             }

        }

        typedef asctoolkit::slam::Map<int> TriCountMap;
        typedef asctoolkit::slam::Map<int> CardinalityVTMap;

        TriCountMap triCount( &m_elementSet);
        for(int lev=0; lev< this->m_levels.size(); ++lev)
        {
            const LeavesLevelMap& levelLeafMap = this->m_leavesLevelMap[ lev ];
            levelBlocks[lev] = levelLeafMap.size();
            for(LeavesIterator it=levelLeafMap.begin(), itEnd = levelLeafMap.end(); it != itEnd; ++it)
            {
                const InOutBlockData& leafData = it->second;
                if(leafData.isLeaf() && leafData.hasData())
                {
                    BlockIndex blk(it->first, lev);
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

        typedef quest::BoundingBox<double,1> MinMaxRange;
        typedef MinMaxRange::PointType LengthType;

        typedef std::map<int,MinMaxRange> LogRangeMap;
        LogRangeMap triCountRange;



        for ( int i=0; i < m_elementSet.size(); ++i )
        {
            LengthType count( triCount[i] );
            int expBase2;
            std::frexp( triCount[i], &expBase2);
            triCountHist[ expBase2 ]++;
            triCountRange[ expBase2 ].addPoint( count );


//            if( triCount[i] > 180)
//                SLIC_INFO( "-- Triangle " << i << " had a bucket count of " << triCount[i] );
        }

        std::stringstream triCountStr;
        triCountStr<<"\tTriangle index count (lg-arithmic): ";
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

        CardinalityVTMap cardVT(&m_vertexSet);
        for ( int i=0; i < m_elementSet.size(); ++i )
        {
            TriVertIndices tvRel = triangleVertexIndices(i);
            cardVT[tvRel[0]]++;
            cardVT[tvRel[1]]++;
            cardVT[tvRel[2]]++;
        }

        typedef std::map<int,int> LinHistogram;
        LinHistogram vtCardHist;
        for ( int i=0; i < m_vertexSet.size(); ++i )
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




//
//        // Add field to the triangle mesh
//        meshtk::FieldData* CD = m_surfaceMesh->getCellFieldData();
//        CD->addField( new meshtk::FieldVariable< TriangleIndex >("blockCount", meshTris.size()) );
//
//        int* blockCount = CD->getField( "blockCount" )->getIntPtr();
//
//        SLIC_ASSERT( blockCount != ATK_NULLPTR );
//
//        for ( int i=0; i < meshTris.size(); ++i ) {
//            blockCount[i] = triCount[i];
//        }
//
//        // Add field to the triangle mesh
//        meshtk::FieldData* ND = m_surfaceMesh->getNodeFieldData();
//        ND->addField( new meshtk::FieldVariable< int >("vtCount", m_vertexSet.size()) );
//
//        int* vtCount = ND->getField( "vtCount" )->getIntPtr();
//
//        SLIC_ASSERT( vtCount != ATK_NULLPTR );
//
//        for ( int i=0; i < m_vertexSet.size(); ++i ) {
//            vtCount[i] = cardVT[i];
//        }
    }
}

template<int DIM>
void InOutOctree<DIM>::checkValid() const
{
#ifdef ATK_DEBUG
    typedef typename OctreeBaseType::MapType LeavesLevelMap;
    typedef typename LeavesLevelMap::const_iterator LeavesIterator;


    std::string treeTypeStr = (m_generationState >= INOUTOCTREE_ELEMENTS_INSERTED)
                ? "PM" : "PR";

    SLIC_DEBUG("Inside InOutOctree::checkValid() to verify state of " << treeTypeStr << " octree.");


    // Iterate through the tree
    // Internal blocks should not have associated vertex data
    // Leaf block consistency depends on 'color'
    //      Black or White blocks should have no vertex data (NO_VERTEX) and no triangles
    //      Gray blocks should have a vertex reference; it may or may not be located within the block
    //          The sum of vertices located within a block should equal the number of mesh vertices.
    //      Gray blocks should have one or more triangles.
    //          All triangles should be incident in a common vertex -- which equals the indexed vertex, if it exists.

    SLIC_ASSERT_MSG( !m_vertexSet.empty()
                   , "InOutOctree::checkValid assumes that we have already inserted "
                     << "the vertices into the octree.");




    if(m_generationState > INOUTOCTREE_VERTICES_INSERTED )
    {
      SLIC_DEBUG("--Checking that each vertex is in a leaf block of the tree.");
      const int numVertices = m_vertexSet.size();
      for(int i=0; i< numVertices; ++i)
      {
          // Check that we can find the leaf block indexing each vertex from its position
            BlockIndex vertBlock = this->findLeafBlock( vertexPosition(i) );
            const InOutBlockData& leafData = (*this)[vertBlock];
            VertexIndex vertInBlock = leafVertex(vertBlock, leafData);

            SLIC_ASSERT_MSG( leafData.hasData() &&  vertInBlock == i
                         ,  " Vertex " << i << " at position " << vertexPosition(i)
                         << " \n\t was not indexed in block " << vertBlock
                         << " with bounding box " << this->blockBoundingBox(vertBlock)
                         << " ( point is" << (this->blockBoundingBox(vertBlock).contains(vertexPosition(i))? "" : " NOT" )
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
                       << "\n\t -- vertex set size: " << m_vertexSet.size()
                       << "-- set is empty " << ( m_vertexSet.empty() ? "yes" : "no")
                        );
      }
    }

    if(m_generationState >= INOUTOCTREE_ELEMENTS_INSERTED)
    {
        SLIC_DEBUG("--Checking that each triangle is referenced by the leaf blocks containing its vertices.");

        const int numTriangles = m_elementSet.size();
        for(int tIdx=0; tIdx< numTriangles; ++tIdx)
        {
            TriVertIndices tvRel = triangleVertexIndices( tIdx );
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
            const LeavesLevelMap& levelLeafMap = this->m_leavesLevelMap[ lev ];
            for(LeavesIterator it=levelLeafMap.begin(), itEnd = levelLeafMap.end(); it != itEnd; ++it)
            {
                const BlockIndex block(it->first, lev);
                const InOutBlockData& data = it->second;

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

                            // Check that the blocks vertex is one of this triangle's vertices
                            TriVertIndices tvRel = triangleVertexIndices( tIdx );
                            SLIC_ASSERT_MSG( incidentInVertex(tvRel, vIdx)
                                             , "All triangles in a gray block must be incident on a common vertex,"
                                             << " but triangles " << tIdx << " with vertices " << tvRel
                                             << " in block " << block << " is not incident in vertex " << vIdx);

                            // Check that this triangle intersects the bounding box of the block
                            GeometricBoundingBox blockBB = this->blockBoundingBox(block);
                            SLIC_ASSERT_MSG( blockIndexesVertex(tIdx, block)
                                             || intersect( trianglePositions(tIdx), blockBB)
                                           , "Triangle " << tIdx << " was indexed in block " << block
                                            << " but it does not intersect the block."
                                            << "\n\tBlock bounding box: " << blockBB
                                            << "\n\tTriangle positions: " << trianglePositions(tIdx)
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
            const LeavesLevelMap& levelLeafMap = this->m_leavesLevelMap[ lev ];
            for(LeavesIterator it=levelLeafMap.begin(), itEnd = levelLeafMap.end(); it != itEnd; ++it)
            {
                const BlockIndex block(it->first, lev);
                const InOutBlockData& data = it->second;

                if(data.isLeaf())
                {
                    SLIC_ASSERT_MSG( data.color() != InOutBlockData::Undetermined
                                , "Problem in block " << block
                                 << " with data " << data
                                 <<" -- Block is uncolored.");

                }
            }
        }
    }

#endif

    SLIC_DEBUG("done.");

}



} // end namespace quest

#endif  // SPATIAL_OCTREE__HXX_
