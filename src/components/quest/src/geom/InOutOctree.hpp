#ifndef INOUT_OCTREE__HXX_
#define INOUT_OCTREE__HXX_

#include "quest/BoundingBox.hpp"
#include "quest/Intersection.hpp"
#include "quest/Point.hpp"
#include "quest/Vector.hpp"
#include "quest/SquaredDistance.hpp"
#include "quest/fuzzy_compare.hpp"


#include "slic/slic.hpp"

#include "slam/Map.hpp"
#include "slam/RangeSet.hpp"


#include "quest/Mesh.hpp"
#include "quest/SpatialOctree.hpp"
#include "quest/Triangle.hpp"
#include "quest/Mesh.hpp"
#include "quest/UnstructuredMesh.hpp"


#include <vector> // For LeafData triangle lists -- TODO replace with SLAM DynamicVariableRelation...
#include <ostream>
#include <set>


#ifdef USE_CXX11
  // Note: Replace the explicit timer functionality once
  //       a timer class is added to common.
  #include <chrono>
  #include <ctime>
#endif


#define DEBUG_VERT_IDX  -2
#define DEBUG_TRI_IDX   -2

#ifndef DUMP_VTK_MESH
//    #define DUMP_VTK_MESH
#endif

/**
 * \file
 *
 * \brief Defines an InOutOctree for containment queries on a surface.
 */


namespace quest
{

    /**
     * \brief DataType for InOut octree leaf nodes
     * Contains a reference to a vertex ID in the mesh
     * (and has a unique ID from its parent class)
     */
    class InOutLeafData : public LeafData
    {
    public:
        enum { NO_VERTEX = -1, BLOCK_COLOR_WHITE = -2, BLOCK_COLOR_BLACK = -3 };
        enum LeafColor { White = -1, Gray = 0, Black = 1};

        typedef int VertexIndex;
        typedef int TriangleIndex;

        typedef std::vector<TriangleIndex> TriangleList;

    public:
        /**
         * \brief Constructor for an InOutOctree LeafData
         * \param vInd The index of a vertex (optional; default is to not set a vertex)
         */
        InOutLeafData(VertexIndex vInd = NO_VERTEX)
            : LeafData()
            , m_vertIndex(vInd)
        {}

        /**
         * \brief Copy constructor for an InOutLeafData instance
         */
        InOutLeafData(const InOutLeafData& other)
            : LeafData(other)
            , m_vertIndex(other.m_vertIndex)
            , m_tris(other.m_tris)
        {
        }

        /**
         * \brief Assignment operator for an InOutLeafData instance
         */
        InOutLeafData& operator=(const InOutLeafData& other)
        {
            LeafData::operator=(other);
            this->m_vertIndex = other.m_vertIndex;
            this->m_tris = other.m_tris;

            return *this;
        }

        /**
         * \brief Checks whether there is a vertex associated with this leaf
         */
        bool hasVertex() const { return m_vertIndex >= 0; }

        /**
         * \brief Check whether this Leaf has any associated triangles
         */
        bool hasTriangles() const { return !m_tris.empty(); }

        /**
         * \brief Find the 'color' of this LeafBlock
         * Black indicates that the entire block is within the surface
         * White indicates that the entire block is outside the surface
         * Gray indicates that the block intersects the surface geometry
         */
        LeafColor color() const
        {
            if(hasVertex() )
            {
                return Gray;
            }
            else if( m_vertIndex == BLOCK_COLOR_BLACK)
            {
                return Black;
            }
            else // WHITE
            {
                SLIC_ASSERT(m_vertIndex == BLOCK_COLOR_WHITE);

                return White;
            }
        }

        /**
         * \brief Find the number of triangles associated with this leaf
         */
        int numTriangles() const { return m_tris.size(); }

        /**
         * \brief Removes all indexed data from this leaf
         */
        void clear()
        {
            m_vertIndex = NO_VERTEX;
            m_tris.clear();
        }

        /**
         * \brief Sets the vertex associated with this leaf
         */
        void setVertex(VertexIndex vInd) { m_vertIndex = vInd; }

        void addTriangle(TriangleIndex tInd) { m_tris.push_back(tInd); }


        /**
         * \brief Accessor for the index of the vertex associated with this leaf
         */
        VertexIndex& vertexIndex() { return m_vertIndex; }

        /**
         * \brief Constant accessor for the index of the vertex associated with this leaf
         */
        const VertexIndex& vertexIndex() const { return m_vertIndex; }


        const TriangleList& triangles() const { return m_tris; }
        TriangleList& triangles() { return m_tris;}


        /**
         * \brief Equality operator to determine if two InOutLeafData instances are equivalent
         */
        friend bool operator==(const InOutLeafData& lhs, const InOutLeafData& rhs )
        {
            return (static_cast<const LeafData&>(lhs) == static_cast<const LeafData&>(rhs))
                && (lhs.m_vertIndex == rhs.m_vertIndex)
                && (lhs.m_tris.size() == lhs.m_tris.size())     // Note: We are not checking the contents
                // && (lhs.m_tris == lhs.m_tris)                //       of the triangle array, only the size
                ;
        }

    private:
        VertexIndex m_vertIndex;
        TriangleList m_tris;
    };

/**
 * \class
 * \brief Handles generation of the spatial index on a surface mesh
 * for containment queries -- given an arbitrary point in space,
 * is it inside or outside of the surface
 */
template<int DIM>
class InOutOctree : public SpatialOctree<DIM, InOutLeafData>
{
public:

    typedef OctreeBase<DIM, InOutLeafData> OctreeBaseType;
    typedef SpatialOctree<DIM, InOutLeafData> SpatialOctreeType;

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
    typedef Point<VertexIndex,NUM_TRI_VERTS> TriVertIndices;

    typedef asctoolkit::slam::PositionSet MeshVertexSet;
    typedef asctoolkit::slam::Map<BlockIndex> VertexBlockMap;


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
     * \brief Helper function to insert a triangle into the octree
     * \param idx The index of the triangle that we are inserting
     * \param triBB A spatial bounding box containing the triangle
     * \param block The octree BlockIndex to which we are considering adding the triangle
     */
    void insertTriangle(TriangleIndex idx
                      , const GeometricBoundingBox& triBB
                      , const BlockIndex& block);


    /**
     * Use octree index over mesh vertices to convert the 'triangle soup'
     * from the stl file into an indexed triangle mesh representation.
     * In particular, all vertices in the mesh that are nearly coincident will be merged,
     * and degenerate triangles (where the three vertices do not have unique indices)
     * will be removed.
     */
    void updateSurfaceMeshVertices();

    /**
     * \brief Helper function to retrieve the position of the vertex from the mesh
     * \param idx The index of the vertex within the surface mesh
     */
    SpacePt vertexPosition(VertexIndex idx) const
    {
        SpacePt pt;
        m_surfaceMesh->getMeshNode(idx, pt.data() );
        return pt;
    }


    /**
     * Utility function to get the indices of the boundary vertices of a triangle
     * \param idx The index of a triangle within the surface mesh
     */
    TriVertIndices triangleVertexIndices(TriangleIndex idx) const
    {
        TriVertIndices vertIds;
        m_surfaceMesh->getMeshCell(idx, vertIds.data());

        return vertIds;
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
        for(int i=0; i< NUM_TRI_VERTS; ++i)
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

    /**
     * Utility function to print some statistics about the InOutOctree instance
     */
    void printOctreeStats(bool trianglesAlreadyInserted) const;

private:

    void dumpOctreeMeshVTK(const std::string& name) const
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
                const InOutLeafData& leafData = (*this)[block];
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

                leafVertID[leafCount] = vIdx;
                leafVertID_unique[leafCount] = (this->findLeafBlock(pt) == block)? vIdx : InOutLeafData::NO_VERTEX;
                leafTriCount[leafCount] = leafData.numTriangles();
                leafCount++;
            }
        }

        // Add the fields
        meshtk::FieldData* CD = debugMesh->getCellFieldData();
        CD->addField( new meshtk::FieldVariable< VertexIndex >("vertID", leafSet.size()) );
        CD->addField( new meshtk::FieldVariable< VertexIndex >("uniqVertID", leafSet.size()) );
        CD->addField( new meshtk::FieldVariable< int >("triCount", leafSet.size()) );

        VertexIndex* vertID = CD->getField( "vertID" )->getIntPtr();
        VertexIndex* uniqVertID = CD->getField( "uniqVertID" )->getIntPtr();
        int* triCount = CD->getField( "triCount" )->getIntPtr();

        SLIC_ASSERT( vertID != ATK_NULLPTR );
        SLIC_ASSERT( uniqVertID != ATK_NULLPTR );
        SLIC_ASSERT( triCount != ATK_NULLPTR );

        for ( int i=0; i < leafSet.size(); ++i ) {
            vertID[i] = leafVertID[i];
            uniqVertID[i] = leafVertID_unique[i];
            triCount[i] = leafTriCount[i];
        }

        debugMesh->toVtkFile(name);

        delete debugMesh;
        debugMesh = ATK_NULLPTR;
      #endif
    }

    /**
     * \brief Utility function to dump a single element vtk file
     */
    void dumpMeshVTK(const std::string& name, int idx, const BlockIndex& block, const GeometricBoundingBox& blockBB, bool isTri) const
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
      #endif
    }

    /**
     * Deal with numerical degeneracies on triangle boundaries by explicitly checking for
     * block containing triangle vertices
     */
    bool blockIndexesVertex(TriangleIndex idx, const BlockIndex& block) const
    {
        TriVertIndices tVerts = triangleVertexIndices(idx);

        if(!m_vertexSet.empty())
        {
            for(int i=0; i<3; ++i)
            {
                if( block == m_vertexToBlockMap[ tVerts[i]])
                    return true;
            }
        }
        else
        {
            for(int i=0; i<3; ++i)
            {
                if( block == this->findLeafBlock( vertexPosition(tVerts[i])))
                    return true;
            }
        }

        return false;
    }

public:
    /**
     * \brief Checks if all indexed triangles in the block share a common vertex
     * \param leafBlock [in] The current octree block
     * \param leafData [inout] The data associated with this block
     * \note A side effect of this function is that we set the leafData's vertex to the common
     * vertex if one is found
     * \return True, if all triangles indexed by this leaf share a common vertex, false otherwise.
     */
    bool allTrianglesIncidentInCommonVertex(const BlockIndex& leafBlock, InOutLeafData& leafData) const;

    bool expensiveBoundingBoxTriTest() const { return false; }
    bool determineInsideOutsideForGrayCell( const SpacePt & pt, const LeafData&  data) { return false;}

private:
  DISABLE_COPY_AND_ASSIGNMENT(InOutOctree)

    /**
     * \brief Checks internal consistency of the octree representation
     */
    void checkValid(bool trianglesAlreadyInserted) const;

protected:
    SurfaceMesh*& m_surfaceMesh; /** pointer to pointer to allow changing the mesh */

    MeshVertexSet m_vertexSet;
    VertexBlockMap m_vertexToBlockMap;
};



template<int DIM>
void InOutOctree<DIM>::generateIndex ()
{
    // Loop through mesh vertices
    SLIC_INFO("Generating index on mesh vertices.");

#ifdef USE_CXX11
    std::chrono::time_point<std::chrono::system_clock> vStart, vEnd;
    vStart = std::chrono::system_clock::now();
#endif

    // STEP 1 -- Add mesh vertices to octree
    int numMeshVerts = m_surfaceMesh->getMeshNumberOfNodes();
    for(int idx=0; idx < numMeshVerts; ++idx)
    {
        insertVertex(idx);
    }

#ifdef USE_CXX11
    vEnd = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = vEnd-vStart;
    SLIC_INFO("\tInserting vertices took " << elapsed_seconds.count() << " seconds.");

    // start timer for updating mesh
    vStart = std::chrono::system_clock::now();
#endif

    updateSurfaceMeshVertices();

#ifdef USE_CXX11
    vEnd = std::chrono::system_clock::now();
    elapsed_seconds = vEnd-vStart;
    SLIC_INFO("\tUpdating mesh took " << elapsed_seconds.count() << " seconds.");

#endif
    SLIC_INFO("** Octree stats after inserting vertices");

    dumpOctreeMeshVTK("prOctree.vtk");
    printOctreeStats(false);
    checkValid(false);

#ifdef USE_CXX11
    // start timer for triangle insertion
    vStart = std::chrono::system_clock::now();
#endif

    // STEP 2 -- Add mesh triangles to octree
    int numMeshTris = m_surfaceMesh->getMeshNumberOfCells();
    for(int idx=0; idx < numMeshTris; ++idx)
    {
        GeometricBoundingBox triBB = triangleBoundingBox(idx);

#ifndef USE_THIS_ONE
        BlockIndex enclosingBlock =
                // TODO: Add code to find smallest enclosing BlockIndex for triBB
                // For now, we can use the root block of the octree
                this->root();

        insertTriangle(idx, triBB, enclosingBlock );

#else
        int lev = this->containingLevel( triBB.range());
        BlockIndex enclosingBlock(this->findGridCellAtLevel(triBB.getMin(), lev), lev);

        BoundingBox<int, 3> gridBB;
        gridBB.addPoint( this->findGridCellAtLevel(triBB.getMin(), lev) );
        gridBB.addPoint( this->findGridCellAtLevel(triBB.getMax(), lev) );
        TriVertIndices tv = triangleVertexIndices(idx);
        for(int i=0; i<3; ++i)
        {
            BlockIndex vBlock = m_vertexToBlockMap[ tv[i]];
            while( vBlock.level() > lev )
                vBlock = vBlock.parent();

            gridBB.addPoint(vBlock.pt());
        }

        const GridPt& minPt = gridBB.getMin();
        const Vector<int,3> range = gridBB.range();

        if(idx < 10)
            SLIC_INFO("-- Triangle " << idx << " with BB -- " << triBB <<" -- covers grid cells -- " << gridBB);

        std::set<BlockIndex> parentBlocks;
        for(int i=0; i<= range[0]; ++i)
            for(int j=0; j<= range[1]; ++j)
                for(int k=0; k<= range[2]; ++k)
                {
                    BlockIndex blk(GridPt::make_point(minPt[0]+i, minPt[1]+j,  minPt[2]+k), lev);
                    while(!this->hasBlock(blk))
                        blk = blk.parent();
                    if(blk.level())
                    {
                        if(parentBlocks.find(blk) != parentBlocks.end())
                            continue;
                        else
                            parentBlocks.insert(blk);
                    }

                    insertTriangle(idx, triBB, blk);
                }
#endif

    }

#ifdef USE_CXX11
    vEnd = std::chrono::system_clock::now();
    elapsed_seconds = vEnd-vStart;
    SLIC_INFO("\tInserting triangles took " << elapsed_seconds.count() << " seconds.");
#endif

    dumpOctreeMeshVTK("pmOctree.vtk");

    SLIC_INFO("** Octree stats after inserting triangles");
    printOctreeStats(true);
    checkValid(true);

    // STEP 3 -- Color the blocks of the octree -- Black (in), White(out), Gray(Intersects surface)


}


template<int DIM>
void InOutOctree<DIM>::insertVertex (VertexIndex idx, int startingLevel)
{
    static const double EPS = 1e-18;


    SpacePt pt = vertexPosition(idx);

//    static GeometricBoundingBox debugBBox = (SpacePt::make_point(0,-26.6102,42.9546));
//    static bool expanded = false;
//    if(!expanded)
//    {
//        debugBBox.expand(0.01);
//        expanded = true;
//    }
//
//    if( debugBBox.contains(pt))
//    {
//        SLIC_INFO("insertVertex -- for vertex " << idx << " at location " << pt << " contained in bbox " << debugBBox );
//    }


    BlockIndex block = this->findLeafBlock(pt, startingLevel);
    InOutLeafData& leafData = (*this)[block];

    if( idx == DEBUG_VERT_IDX)
    {
        SLIC_INFO("\t -- pt coords: " << pt
                  << " contained in block is " << block
                  << " blockBB " << this->blockBoundingBox(block)
                  << " leaf vertex is " << leafData.vertexIndex()
                    );
    }


    if(! leafData.hasVertex())
    {
        leafData.setVertex(idx);

        // Update the vertex to block map for this vertex
        if(!m_vertexSet.empty())
            m_vertexToBlockMap[idx] = block;

    }
    else
    {
        // check if we should merge the vertices
        VertexIndex origVertInd = leafData.vertexIndex();
        if( squared_distance( pt, vertexPosition(origVertInd) ) >= EPS )
        {
            leafData.clear();
            this->refineLeaf(block);

            insertVertex(origVertInd, block.childLevel() );
            insertVertex(idx, block.childLevel() );
        }

    }

//    if( debugBBox.contains(pt) && leafData.vertexIndex() == idx)
    if( leafData.vertexIndex() == DEBUG_VERT_IDX)
    {
        SLIC_INFO("\t -- vertex " << idx << " is indexed in block " << block);
    }


}


template<int DIM>
void InOutOctree<DIM>::updateSurfaceMeshVertices()
{
    typedef asctoolkit::slam::PositionSet MeshVertsSet;
    typedef asctoolkit::slam::Map<VertexIndex> IndexMap;

    // Create a map from old vertex IDs to new vertex ids
    MeshVertsSet oldVerts( m_surfaceMesh->getMeshNumberOfNodes() );
    IndexMap vertexIndexMap( &oldVerts, InOutLeafData::NO_VERTEX);

    // Generate unique indices for mesh vertices
    int uniqueVertexCounter = 0;
    for(int i=0; i< oldVerts.size(); ++i)
    {
        // Find the block and its indexed vertex in the octree
        BlockIndex leafBlock = this->findLeafBlock( vertexPosition(i) );
        SLIC_ASSERT( (*this)[leafBlock].hasVertex() );
        VertexIndex vInd = (*this)[ leafBlock ].vertexIndex();

        // If the indexed vertex doesn't have a new id, give it one
        if(vertexIndexMap[vInd] == InOutLeafData::NO_VERTEX)
        {
            vertexIndexMap[vInd] = uniqueVertexCounter++;
        }

        // If this is not the indexed vertex, grab that vertex's new IDX
        if(vInd != i)
        {
            vertexIndexMap[i] = vertexIndexMap[vInd];
        }
    }

    // Find coordinates for the new mesh vertices
    // Cache in temporary SLAM Map since these may be out of order w.r.t. new vertex indices
    typedef asctoolkit::slam::Map<SpacePt> MeshCoordsField;

    m_vertexSet = MeshVertexSet( uniqueVertexCounter );

    MeshCoordsField newVertCoords( &m_vertexSet);
    for(int i=0; i< oldVerts.size(); ++i)
    {
        const VertexIndex vInd = vertexIndexMap[i];
        newVertCoords[vInd] = vertexPosition(i);
    }



    // Create a new triangles mesh.
    // * add new vertices from octree to the mesh
    // * update vertex references in octree to new vertex indices
    typedef meshtk::UnstructuredMesh< meshtk::LINEAR_TRIANGLE > TriangleMesh;
    TriangleMesh* newMesh = new TriangleMesh(3);
    m_vertexToBlockMap = VertexBlockMap(&m_vertexSet);
    for(int i = 0; i< m_vertexSet.size(); ++i)
    {
        const SpacePt& pos = newVertCoords[i];
        newMesh->insertNode(pos[0], pos[1], pos[2]);

        BlockIndex leafBlock = this->findLeafBlock(pos);
        SLIC_ASSERT( this->isLeaf(leafBlock) && (*this)[leafBlock].hasVertex() );
        (*this)[ leafBlock ].setVertex(i);
        m_vertexToBlockMap[i] = leafBlock;
    }

    // Add triangles from old mesh to new mesh using updated vertex ids
    int numOldMeshTris =  m_surfaceMesh->getMeshNumberOfCells();
    for(int i=0; i< numOldMeshTris ; ++i)
    {
       TriVertIndices vertIndices = triangleVertexIndices(i);

       // Remap the vertex IDs
       for(int j=0; j< 3; ++j)
           vertIndices[j] = vertexIndexMap[ vertIndices[j] ];

       // Skip degenerate triangles -- need 3 unique vertex IDS
       if(    (vertIndices[0] != vertIndices[1])
           && (vertIndices[1] != vertIndices[2])
           && (vertIndices[2] != vertIndices[0]) )
       {
           newMesh->insertCell( vertIndices.data(), meshtk::LINEAR_TRIANGLE, 3);
       }
    }

    // Delete old mesh, redirect pointer to newly created mesh
    delete m_surfaceMesh;
    m_surfaceMesh = newMesh;
}


template<int DIM>
void InOutOctree<DIM>::insertTriangle (TriangleIndex idx
                                     , const GeometricBoundingBox& triBB
                                     , const BlockIndex& block)
{
    GeometricBoundingBox blockBB = this->blockBoundingBox(block);

    std::string indent;

    if(idx == DEBUG_TRI_IDX)
    {
        std::stringstream sstr;
        for(int i=0; i< block.level(); ++i)
            sstr<<"  ";
        indent = sstr.str();
    }

    if(idx == DEBUG_TRI_IDX)
    {
        static bool triDumped = false;
        if(!triDumped)
        {
            // Dump the triangle
            dumpMeshVTK("triangle", idx, block, blockBB, true);

            // Dump the bounding boxes of its three vertices
            TriVertIndices vertInd = triangleVertexIndices(idx);
            for(int i=0;i<3;++i)
            {
                SpacePt pt = vertexPosition( vertInd[i]);
                BlockIndex blk = this->findLeafBlock(pt );
                GeometricBoundingBox blkBB = this->blockBoundingBox( blk);

                SLIC_INFO(" ** Vertex " << i << " of triangle " << idx
                        << " -- at position " << pt
                        << " -- located in block " << blk
                        << " -- with bb " << blkBB
                        );

                dumpMeshVTK("triangleVertBox", idx, blk, blkBB, false);
            }
            triDumped = true;

            bool hasMissingTris = false;
            if(hasMissingTris)
            {
                // Dump any extra problematic bounding box
                std::vector<BlockIndex> missingBlocks;
                missingBlocks.push_back( BlockIndex( Point<int, DIM>::make_point(16,31,16), 5) );
                missingBlocks.push_back( BlockIndex( Point<int, DIM>::make_point(8,15,8), 4) );

                for(size_t i=0; i< missingBlocks.size(); ++i)
                {
                    const BlockIndex& blk = missingBlocks[i];
                    GeometricBoundingBox blkBB = this->blockBoundingBox( blk);

                    SLIC_INFO( " -- block " << blk
                               << " with BB " << blkBB
                               << ( intersect(trianglePositions(idx), blkBB) ? " " : "does not " ) << " intersect triangle " << idx);
                    dumpMeshVTK("triangleMissing", idx, blk, blkBB, false);
                }
            }
        }
    }



    if(idx == DEBUG_TRI_IDX )
    {
        BlockIndex blk( Point<int,DIM>::make_point(16,31,16), 5);
        if( block == blk)
        {
            SLIC_INFO(  "--> Checking triangle " << idx << " against block " << block
                        << "; triangle BB is: " << triBB << "; blockBB is: " << blockBB);
        }
        else if(block == blk.parent())
        {
            SLIC_INFO(  "--> Checking triangle " << idx << " against block " << block
                        << " (parent of blk == " << blk << ")"
                        << "; triangle BB is: " << triBB << "; blockBB is: " << blockBB);
        }
    }


    // Simple initial check -- does triangles' BB overlap with block's BB
    if( intersect( triBB, blockBB) )
    {
//        if( idx == DEBUG_TRI_IDX)
//        {
//            SLIC_INFO("Insert vertex for triangle " << idx
//                      << " against block " << block
//                      << "\n\t\t -- tri -- " << trianglePositions(idx)
//            );
//        }

        if( this->isInternal(block) )
        {
            // Recursive call to children for internal blocks
            for(int i=0; i< block.numChildren(); ++i)
            {
                if(this->hasBlock(block))
                    insertTriangle(idx, triBB, block.child(i));
            }
        }
        else
        {
            SLIC_ASSERT( this->hasBlock(block) && this->isLeaf(block));

            if(idx == DEBUG_TRI_IDX)
            {
                SLIC_INFO( indent << "--> Checking triangle " << idx << " against leaf block " << block
                          << "; triangle BB is: " << triBB << "; blockBB is: " << blockBB);
            }

          if( blockIndexesVertex(idx, block) || intersect( trianglePositions(idx), blockBB) ) // internal block -- perform expensive test
          {
            if(idx == DEBUG_TRI_IDX)
              dumpMeshVTK("triangleIntersect", idx, block, blockBB, false);

            InOutLeafData& leafData = (*this)[block];
            leafData.addTriangle(idx);

            if(idx == DEBUG_TRI_IDX)
                SLIC_INFO(indent << "\t\t-- triangle " << idx << " intersects leaf " << block << " -- added to leafData with " << leafData.numTriangles() << " triangles.");


            // Here we have our refinement criterion for triangles
            // All indexed triangles must be incident in a single common vertex
            if( ! allTrianglesIncidentInCommonVertex(block, leafData) )//&& false)
            {
                // Grab the triangle data
                InOutLeafData::TriangleList parentTriangles;
                parentTriangles.swap(leafData.triangles());
                int numTriangles = parentTriangles.size();

                VertexIndex vIdx = leafData.vertexIndex();

                if( vIdx == DEBUG_VERT_IDX && DEBUG_VERT_IDX >= 0)
                {
                    SLIC_INFO("** insertTriangle -- vertex " << vIdx
                              << " @ position -- " << vertexPosition(vIdx)
                              << " -- parent block " << block
                              << " -- blockBB " << blockBB
                              << " -- vertex is " << ( blockBB.halfOpenContains( vertexPosition(vIdx))? "":"not ")
                              << "contained in block."
                              );

                    for(int i=0; i< 8; ++i)
                    {
                        GeometricBoundingBox childBB = this->blockBoundingBox( block.child(i) );
                        SLIC_INFO("\t-- vertex " << vIdx << " is " << ( childBB.halfOpenContains( vertexPosition(vIdx))? "":"not ")
                                  << " contained in child " << i << " with block " << block.child(i)
                                  << " with bb " <<  childBB
                                  );
                    }
                }


                // Clear the current leaf and refine
                leafData.clear();
                this->refineLeaf(block);

                if( m_vertexToBlockMap[vIdx] == block )
                    insertVertex(vIdx, block.childLevel());

                // Reinsert triangles in children
                for(int i=0; i< numTriangles; ++i)
                {
                    ///  NOTE: Recomputing bounding box might be too expensive --
                    TriangleIndex tIdx = parentTriangles[i];
                    GeometricBoundingBox tBB = triangleBoundingBox(tIdx);

                    if(tIdx == DEBUG_TRI_IDX)
                        SLIC_INFO(indent << "\t\t\t-- removing triangle " << tIdx << " from refining leaf " << block );

                    for(int j=0; j< block.numChildren(); ++j)
                    {
                        insertTriangle(tIdx, tBB, block.child(j));
                    }
                }

//                if(idx == DEBUG_TRI_IDX)
//                    SLIC_INFO(indent << "\t\t-- end insert triangle " << idx << " into leaf " << block );

            }
          }

          else {
                  if(idx == DEBUG_TRI_IDX)
                  {
                      SLIC_INFO( indent << "xx> triangle " << idx << " did not intersect " << block
                                << "; triangle BB is: " << triBB << "; blockBB is: " << blockBB);
                  }
          }
        }
    }
}


template<int DIM>
bool InOutOctree<DIM>::allTrianglesIncidentInCommonVertex(const BlockIndex& leafBlock
                                                         , InOutLeafData& leafData) const
{
    bool shareCommonVert = false;
    VertexIndex commonVert = InOutLeafData::NO_VERTEX;

    const int numTris = leafData.numTriangles();

    if(numTris == 0)
    {
        // No need to refine here
        shareCommonVert = true;
    }
    else if(leafData.hasVertex() && m_vertexToBlockMap[ leafData.vertexIndex()] == leafBlock )
    {

        // Loop through the triangles and check that all are incident against this vertex
        commonVert = leafData.vertexIndex();;

        const InOutLeafData::TriangleList& tris = leafData.triangles();


//        if( commonVert == DEBUG_VERT_IDX)
//        {
//            SLIC_INFO("allTrianglesCommon -- vertex " << commonVert);
//        }


        bool hasTriMissingVert = false;
        for(int i=0; !hasTriMissingVert && i< numTris; ++i)
        {
            TriVertIndices tvRel = triangleVertexIndices(tris[i]);
            if( ! ( (tvRel[0] == commonVert) || (tvRel[1] == commonVert) || (tvRel[2] == commonVert) ) )
            {
                hasTriMissingVert = true;
            }
        }

        if(!hasTriMissingVert)
        {
            shareCommonVert = true;
        }
    }
    else
    {
        // Simple strategy -- can make more efficient later
        // Note: Assumption is that we have very few (i.e. 6 or fewer) triangles in a bucket
        //       Optimize for cases where there are only 1 or 2 indexed triangles
        if(numTris == 1)
        {
            commonVert = triangleVertexIndices(leafData.triangles()[0])[0];
            shareCommonVert = true;
        }
        else if(numTris == 2)
        {
            const InOutLeafData::TriangleList& tris = leafData.triangles();
            TriVertIndices tvRel0 = triangleVertexIndices(tris[0]);
            TriVertIndices tvRel1 = triangleVertexIndices(tris[1]);

            for(int i=0; !shareCommonVert && i<3; ++i)
            {
                commonVert = tvRel1[i];

                if( tvRel0[0] == commonVert || tvRel0[1] == commonVert || tvRel0[2] == commonVert)
                    shareCommonVert = true;
            }

        }
        else    // numTris > 3
        {

            //  1. add verts of all but last triangle to histogram
            std::map<VertexIndex, int> vertCount;
            const InOutLeafData::TriangleList& tris = leafData.triangles();

            const int nTrisMinusOne = numTris -1;
            for(int i=0; i< nTrisMinusOne; ++i)
            {
                TriVertIndices tvRel = triangleVertexIndices(tris[i]);
                ++vertCount[tvRel[0]];
                ++vertCount[tvRel[1]];
                ++vertCount[tvRel[2]];
            }

            //  2. Check that at least one vertex in last triangle has the right count
            //      Note: this correctly handles the case where there is only one triangle
            TriVertIndices tvRel = triangleVertexIndices(tris[nTrisMinusOne]);
            for(int i=0; !shareCommonVert && i<3; ++i)
            {
                VertexIndex vIdx = tvRel[i];
                if(vertCount[vIdx] == nTrisMinusOne)
                {
                    commonVert = vIdx;
                    shareCommonVert = true;
                }
            }
        }
    }

    leafData.setVertex( commonVert );
    return shareCommonVert;
}


template<int DIM>
bool InOutOctree<DIM>::within(const SpacePt& pt) const
{
    if(blockBoundingBox( this->root() ).contains(pt) )
    {
        bool isWithin = false;

        BlockIndex block = this->findLeafBlock(pt);
        InOutLeafData& data = (*this)[block];

        switch( data.color() )
        {
        case InOutLeafData::LeafColor::Black:
            isWithin = true;
            break;
        case InOutLeafData::LeafColor::White:
            isWithin  = false;
            break;
        case InOutLeafData::LeafColor::Gray:
            isWithin = determineInsideOutsideForGrayCell( pt, data);
            break;
        }

        return isWithin;
    }

    return false;
}


template<int DIM>
void InOutOctree<DIM>::printOctreeStats(bool trianglesAlreadyInserted) const
{
    typedef typename OctreeBaseType::MapType LeavesLevelMap;
    typedef typename LeavesLevelMap::const_iterator LeavesIterator;


    typedef asctoolkit::slam::Map<int> LeafCountMap;

    LeafCountMap levelBlocks( &this->m_levels);
    LeafCountMap levelLeaves( &this->m_levels);
    LeafCountMap levelLeavesWithVert( &this->m_levels);
    LeafCountMap levelTriangleRefCount( &this->m_levels);

    int totalBlocks = 0;
    int totalLeaves = 0;
    int totalLeavesWithVert = 0;
    int totalTriangleRefCount = 0;

    // Iterate through blocks -- count the numbers of internal and leaf blocks
    //
    for(int lev=0; lev< this->m_levels.size(); ++lev)
    {
        const LeavesLevelMap& levelLeafMap = this->m_leavesLevelMap[ lev ];
        levelBlocks[lev] = levelLeafMap.size();
        for(LeavesIterator it=levelLeafMap.begin(), itEnd = levelLeafMap.end(); it != itEnd; ++it)
        {
            if(it->second.isLeaf())
                ++levelLeaves[lev];

            if(it->second.hasVertex())
                ++levelLeavesWithVert[ lev ];

            if(it->second.hasTriangles())
                levelTriangleRefCount[ lev ] += it->second.numTriangles();
        }

        totalBlocks += levelBlocks[lev];
        totalLeaves += levelLeaves[lev];
        totalLeavesWithVert += levelLeavesWithVert[lev];
        totalTriangleRefCount += levelTriangleRefCount[ lev ];
    }


    std::stringstream octreeStatsStr;
    octreeStatsStr << "*** " << (trianglesAlreadyInserted ? "PM" : "PR") << " octree summary *** \n";

    for(int lev=0; lev< this->m_levels.size(); ++lev)
    {
        if(levelBlocks[lev] > 0)
        {
            octreeStatsStr << "\t Level " << lev
                           << " has " << levelLeavesWithVert[lev] << " leaves with vert"
                           << " out of " << levelLeaves[lev] << " leaves; "
                           <<  levelBlocks[lev] - levelLeaves[lev] << " internal blocks."
                           << " and " << levelTriangleRefCount[lev] << " triangle references.\n";
        }
    }


    double meshNumTriangles = static_cast<double>(m_surfaceMesh->getMeshNumberOfCells());
    octreeStatsStr<<"  Mesh has " << m_surfaceMesh->getMeshNumberOfNodes() << " vertices."
                 <<"\n  Octree has " << totalLeavesWithVert << " filled leaves; "
                 <<  totalLeaves << " leaves; "
                 <<  totalBlocks - totalLeaves << " internal blocks;"
                 << totalBlocks << " overall blocks."
                 <<" \n\t There were " << totalTriangleRefCount << " triangle references "
                 <<" (avg. " << totalTriangleRefCount /  meshNumTriangles << " refs per triangle).\n";

    SLIC_INFO( octreeStatsStr.str() );

    if(trianglesAlreadyInserted)
    {
        typedef asctoolkit::slam::Map<int> TriCountMap;
        typedef asctoolkit::slam::Map<int> CardinalityVTMap;

        typedef asctoolkit::slam::PositionSet MeshTriangleSet;
        MeshTriangleSet meshTris( m_surfaceMesh->getMeshNumberOfCells());
        TriCountMap triCount( &meshTris);
        for(int lev=0; lev< this->m_levels.size(); ++lev)
        {
            const LeavesLevelMap& levelLeafMap = this->m_leavesLevelMap[ lev ];
            levelBlocks[lev] = levelLeafMap.size();
            for(LeavesIterator it=levelLeafMap.begin(), itEnd = levelLeafMap.end(); it != itEnd; ++it)
            {
                const InOutLeafData& leafData = it->second;
                if(leafData.isLeaf() && leafData.hasTriangles())
                {
                    const InOutLeafData::TriangleList& tris = leafData.triangles();
                    for(int i = 0; i < leafData.numTriangles(); ++i)
                    {
                        ++triCount[ tris[i]];
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



        for ( int i=0; i < meshTris.size(); ++i )
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
        for ( int i=0; i < meshTris.size(); ++i )
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





        // Add field to the triangle mesh
        meshtk::FieldData* CD = m_surfaceMesh->getCellFieldData();
        CD->addField( new meshtk::FieldVariable< TriangleIndex >("blockCount", meshTris.size()) );

        int* blockCount = CD->getField( "blockCount" )->getIntPtr();

        SLIC_ASSERT( blockCount != ATK_NULLPTR );

        for ( int i=0; i < meshTris.size(); ++i ) {
            blockCount[i] = triCount[i];
        }

        // Add field to the triangle mesh
        meshtk::FieldData* ND = m_surfaceMesh->getNodeFieldData();
        ND->addField( new meshtk::FieldVariable< int >("vtCount", m_vertexSet.size()) );

        int* vtCount = ND->getField( "vtCount" )->getIntPtr();

        SLIC_ASSERT( vtCount != ATK_NULLPTR );

        for ( int i=0; i < m_vertexSet.size(); ++i ) {
            vtCount[i] = cardVT[i];
        }
    }



}



namespace{
    std::ostream& operator<<(std::ostream& os, const std::vector<int>& vec)
    {
        os<<"[";
        for(std::size_t i=0; i< vec.size(); ++i)
            os << vec[i] << ( (i == vec.size()-1) ? "]" : ",");

        return os;
    }
}

template<int DIM>
void InOutOctree<DIM>::checkValid(bool trianglesAlreadyInserted) const
{
#ifdef ATK_DEBUG
    typedef typename OctreeBaseType::MapType LeavesLevelMap;
    typedef typename LeavesLevelMap::const_iterator LeavesIterator;

    SLIC_INFO("Inside InOutOctree::checkValid() to verify state of PR or PM octree.");


    // Iterate through the tree
    // Internal blocks should not have associated vertex data
    // Leaf block consistency depends on color
    //      Black or White blocks should have no vertex data (NO_VERTEX) and no triangles
    //      Gray blocks should have a vertex; it may or may not be located within the block
    //          The sum of vertices within the block should equal the number of mesh vertices.
    //      Gray nodes should have one or more triangles.
    //          All triangles should be incident in a common vertex -- which equals the indexed vertex, if it exists.


    SLIC_INFO("--Checking that each vertex is in a leaf block of the tree.");
    const int numVertices = m_surfaceMesh->getMeshNumberOfNodes();
    for(int i=0; i< numVertices; ++i)
    {
        BlockIndex vertBlock = this->findLeafBlock( vertexPosition(i) );
        const InOutLeafData& leafData = (*this)[vertBlock];
        SLIC_ASSERT_MSG( leafData.hasVertex() && leafData.vertexIndex() == i
                     ,  " Vertex " << i << " at position " << vertexPosition(i)
                     << " \n\t was not indexed in block " << vertBlock
                     << " with bounding box " << this->blockBoundingBox(vertBlock)
                     << " ( point is" << (this->blockBoundingBox(vertBlock).contains(vertexPosition(i))? "" : " NOT" )
                     << " contained in block )."
                     << " \n\n *** \n Leaf data: "
                     << " \n\t Is leaf: " << ( leafData.isLeaf()? " yes" : " no")
                     << " \n\t Has vertex: " << ( leafData.hasVertex()? " yes" : " no")
                     << " \n\t Vertex ID: " << leafData.vertexIndex()
                     << " \n\t Has triangles? " << ( leafData.hasTriangles()? " yes" : " no")
                     << " \n\t Triangle count " << leafData.numTriangles()
                     << " \n\t\t indices: " << leafData.triangles()
                     << " \n ***"
        );

        // Check that our cached value of the vertex's block is accurate
        if(!m_vertexSet.empty())
        {
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

    // Check that each triangle is referenced in at least the blocks indexing its vertices
    if(trianglesAlreadyInserted)
    {
        SLIC_INFO("--Checking that each triangle is referenced by the leaf blocks containing its vertices.");
        const int numTriangles = m_surfaceMesh->getMeshNumberOfCells();
        for(int tIdx=0; tIdx< numTriangles; ++tIdx)
        {
            TriVertIndices tvRel = triangleVertexIndices( tIdx );
            for(int j=0; j<3;++j)
            {
                VertexIndex vIdx = tvRel[j];
                BlockIndex vertBlock = m_vertexToBlockMap[vIdx];
                const InOutLeafData& leafData = (*this)[vertBlock];

                // Check that this triangle is referenced here.
                bool foundTriangle = false;
                const InOutLeafData::TriangleList& leafTris = leafData.triangles();
                for(int k=0; !foundTriangle && k< leafData.numTriangles(); ++k)
                {
                    if( leafTris[k] == tIdx)
                        foundTriangle = true;
                }

                SLIC_ASSERT_MSG( foundTriangle
                               , "Did not find triangle with index " << tIdx << " and vertices" << tvRel
                               << " in block " << vertBlock << " containing vertex " << vIdx
                               << " \n\n *** \n Leaf data: "
                               << " \n\t Is leaf: " << ( leafData.isLeaf()? " yes" : " no")
                               << " \n\t Has vertex: " << ( leafData.hasVertex()? " yes" : " no")
                               << " \n\t Vertex ID: " << leafData.vertexIndex()
                               << " \n\t Has triangles? " << ( leafData.hasTriangles()? " yes" : " no")
                               << " \n\t Triangle count " << leafData.numTriangles()
                               << " \n\t\t indices: " << leafData.triangles()
                               << " \n ***"
                               );
            }
        }

        // Check that internal blocks have no triangle / vertex
        //       and leaf blocks satisfy the conditions above
        SLIC_INFO("--Checking that internal blocks have no data, and that leaves satisfy all PM conditions");
        for(int lev=0; lev< this->m_levels.size(); ++lev)
        {
            const LeavesLevelMap& levelLeafMap = this->m_leavesLevelMap[ lev ];
            for(LeavesIterator it=levelLeafMap.begin(), itEnd = levelLeafMap.end(); it != itEnd; ++it)
            {
                const BlockIndex block(it->first, lev);
                const InOutLeafData& data = it->second;

                if( !data.isLeaf() )
                {
                    SLIC_ASSERT( !data.hasVertex() );
                    SLIC_ASSERT( !data.hasTriangles() );
                }
                else // leaf block
                {
                    if( data.hasVertex())
                    {
                        VertexIndex vIdx = data.vertexIndex();
                        for( int i  = 0; i< data.numTriangles(); ++i)
                        {
                            TriangleIndex tIdx = data.triangles()[i];

                            // Check that the blocks vertex is one of this triangle's vertices
                            TriVertIndices tvRel = triangleVertexIndices( tIdx );
                            SLIC_ASSERT( tvRel[0] == vIdx || tvRel[1] == vIdx || tvRel[2] == vIdx );

                            // Check that this triangle intersects the bounding box of the block
                            GeometricBoundingBox blockBB = this->blockBoundingBox(block);
                            SLIC_ASSERT_MSG( blockIndexesVertex(tIdx, block)
                                             || intersect( trianglePositions(tIdx), blockBB)
                                           , "Triangle " << tIdx << " was indexed in block " << block
                                            << " but it does not intersect the block."
                                            << "\n\tBlock bounding box: " << blockBB
                                            << "\n\tTriangle positions: " << trianglePositions(tIdx)
                                            << "\n\tTriangle vertex indices: " << tvRel
                                            << "\n\tLeaf vertex is: " << data.vertexIndex()
                                            << "\n\tLeaf triangles: " << data.triangles()
                                            << "(" << data.numTriangles() <<")"
                            );
                        }
                    }
                    else
                    {
                        SLIC_ASSERT( !data.hasTriangles() );
                    }
                }

            }
        }
    }


//    static GeometricBoundingBox debugBBox = (SpacePt::make_point(0,-26.6102,42.9546));
//    static bool expanded = false;
//    if(!expanded)
//    {
//        debugBBox.expand(0.01);
//        expanded = true;
//    }
//
//    for(int i=0; i< numVertices; ++i)
//    {
//        SpacePt pt = vertexPosition(i);
//        if( debugBBox.contains(pt) )
//        {
//            BlockIndex vertBlock = this->findLeafBlock( pt );
//            const InOutLeafData& leafData = (*this)[vertBlock];
//
//            SLIC_INFO("*** check -- vertex " << i << " @ coords: " << pt
//                  << " contained in block " << vertBlock
//                  << " w/ blockBB " << this->blockBoundingBox(vertBlock)
//                  << " leaf vertex is " << leafData.vertexIndex()
//                    );
//
//        }
//    }

#endif
}



} // end namespace quest

#endif  // SPATIAL_OCTREE__HXX_
