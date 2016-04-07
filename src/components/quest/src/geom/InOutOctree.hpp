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


#define DEBUG_VERT_IDX  -2
#define DEBUG_TRI_IDX   -2

#ifndef DUMP_VTK_MESH
//    #define DUMP_VTK_MESH
#endif

#ifndef DUMP_OCTREE_INFO
//    #define DUMP_OCTREE_INFO
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
    class InOutLeafData : public BlockData
    {
    public:
        enum { NO_VERTEX = -1, BLOCK_COLOR_WHITE = -2, BLOCK_COLOR_BLACK = -3 };
        enum LeafColor { Undetermined = -2, White = -1, Gray = 0, Black = 1};

        typedef int VertexIndex;
        typedef int TriangleIndex;

        typedef std::vector<TriangleIndex> TriangleList;

    public:
        /**
         * \brief Default constructor for an InOutLeafData
         */
        InOutLeafData()
            : BlockData()
            , m_vertIndex(NO_VERTEX)
        {}

        /**
         * \brief Constructor for an InOutLeafData
         * \param vInd The index of a vertex (optional; default is to not set a vertex)
         */
        InOutLeafData(VertexIndex vInd)
            : BlockData()
            , m_vertIndex(vInd)
        {}

        /**
         * \brief Copy constructor for an InOutLeafData instance
         */
        InOutLeafData(const InOutLeafData& other)
            : BlockData(other)
            , m_vertIndex(other.m_vertIndex)
            , m_tris(other.m_tris)
        {
        }

        /**
         * \brief Assignment operator for an InOutLeafData instance
         */
        InOutLeafData& operator=(const InOutLeafData& other)
        {
            BlockData::operator=(other);
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

        void reserveTriangles(int count) { m_tris.reserve(count); }

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

            switch( m_vertIndex)
            {
            case BLOCK_COLOR_BLACK:     return Black;
            case BLOCK_COLOR_WHITE:     return White;
            default:                    return Undetermined;
            }
        }

        bool isColored() const
        {
            return color() != Undetermined;

        }

        /**
         * \brief Sets the color of the InOutOctree block
         * \pre The block must not already have an associated color
         */
        void setColor(LeafColor leafColor)
        {
            SLIC_ASSERT_MSG( !isColored() || leafColor == color()
                             , "Error setting color of InOutOctreeLeaf."
                             << "\n\t Leaf " << ( isColored() ? " is " : " is not ") << "colored."
                             << "\n\t Passed in leaf color is: " << static_cast<int>(leafColor)
                             << "\n\t Leaf color is " << static_cast<int>( color() )
                             );

            switch(leafColor)
            {
            case Black:  m_vertIndex = BLOCK_COLOR_BLACK; break;
            case White:  m_vertIndex = BLOCK_COLOR_WHITE; break;
            default:
                SLIC_ASSERT_MSG(false, "Error setting color of InOutOctreeLeaf."
                                << " Must be black or white.  Value was "
                                << static_cast<int>(leafColor) );
                break;
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
            m_tris = TriangleList(0);    // reconstruct to deallocate memory
        }

        /**
         * \brief Sets the vertex associated with this leaf
         */
        void setVertex(VertexIndex vInd) { m_vertIndex = vInd; }

        void clearVertex() { m_vertIndex = NO_VERTEX; }

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
            return (static_cast<const BlockData&>(lhs) == static_cast<const BlockData&>(rhs))
                && (lhs.m_vertIndex == rhs.m_vertIndex)
                && (lhs.m_tris.size() == rhs.m_tris.size())     // Note: We are not checking the contents
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

    typedef asctoolkit::slam::PositionSet MeshVertexSet;
    typedef asctoolkit::slam::Map<BlockIndex> VertexBlockMap;

    typedef asctoolkit::slam::PositionSet MeshElementSet;
    typedef asctoolkit::slam::Map<SpacePt> VertexPositionMap;

    typedef asctoolkit::slam::policies::CompileTimeStrideHolder<VertexIndex, NUM_TRI_VERTS>  TVStride;
    typedef asctoolkit::slam::StaticConstantRelation<TVStride, MeshElementSet, MeshVertexSet> TriangleVertexRelation;

    typedef typename TriangleVertexRelation::RelationSet TriVertIndices;


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
     * \brief Insert all triangles of the mesh into the octree, generating a PM octree
     */
    void insertMeshTriangles();


    /**
     * \brief Helper function to insert a vertex into the octree
     * \param idx The index of the vertex that we are inserting
     * \param startingLevel (optional, default = 0) The octree level at which
     * to begin the search for the leaf node covering this vertex
     */
    void insertVertex(VertexIndex idx, int startingLevel = 0);


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

    /**
     * Utility function to print some statistics about the InOutOctree instance
     */
    void printOctreeStats(bool trianglesAlreadyInserted) const;

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

        if(!m_vertexSet.empty())
        {
            for(int i=0; i< tVerts.size(); ++i)
            {
                if( block == m_vertexToBlockMap[ tVerts[i]])
                    return true;
            }
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



    bool colorLeafAndNeighbors(const BlockIndex& blk, InOutLeafData& blkData);


private:
    DISABLE_COPY_AND_ASSIGNMENT(InOutOctree);

    /**
     * \brief Checks internal consistency of the octree representation
     */
    void checkValid(bool trianglesAlreadyInserted) const;


    /**
     * \brief Determines whether the specified point is within the gray leaf
     * \param pt The point we are querying
     * \param data The data associated with the leaf block
     * \return True, if the point is inside the local surface associated with this block, false otherwise
     */
    bool withinGrayBlock( const SpacePt & pt, const InOutLeafData&  data) const;

  /**
    * \brief Add the vertex positions and triangle boundary relations to the surface mesh
    */
    void regenerateSurfaceMesh();

    bool hasGrayChild(const BlockIndex& blk) const;


protected:
    SurfaceMesh*& m_surfaceMesh; /** pointer to pointer to allow changing the mesh */

    MeshVertexSet m_vertexSet;
    VertexBlockMap m_vertexToBlockMap;

    MeshElementSet m_elementSet;
    VertexPositionMap m_vertexPositions;
    TriangleVertexRelation m_triangleToVertexRelation;
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

    /**
     * \brief Utility function to print the triangle indices in a leaf block
     */
    std::ostream& operator<<(std::ostream& os, const std::vector<int>& vec)
    {
        os<<"[";
        for(std::size_t i=0; i< vec.size(); ++i)
            os << vec[i] << ( (i == vec.size()-1) ? "]" : ",");

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
    SLIC_INFO("\tInserting vertices took " << timer.elapsed() << " seconds.");

    // STEP 1(b) -- Update the mesh vertices and cells with after vertex welding from octree
    timer.start();
    updateSurfaceMeshVertices();
    timer.stop();
    SLIC_INFO("\tUpdating mesh took " << timer.elapsed() << " seconds.");


  #ifdef DUMP_OCTREE_INFO
    // -- Print some stats about the octree
    SLIC_INFO("** Octree stats after inserting vertices");
    dumpOctreeMeshVTK("prOctree.vtk", false);
    printOctreeStats(false);
  #endif
    checkValid(false);


    // STEP 2 -- Add mesh triangles to octree
    timer.start();
    insertMeshTriangles();
    timer.stop();

    SLIC_INFO("\tInserting triangles took " << timer.elapsed() << " seconds.");

    // STEP 3 -- Color the blocks of the octree -- Black (in), White(out), Gray(Intersects surface)
    timer.start();
    colorOctreeLeaves();

    timer.stop();
    SLIC_INFO("\tColoring octree leaves took " << timer.elapsed() << " seconds.");

    // -- Print some stats about the octree
  #ifdef DUMP_OCTREE_INFO
    SLIC_INFO("** Octree stats after inserting triangles");
    dumpOctreeMeshVTK("pmOctree.vtk", true);
    printOctreeStats(true);
  #endif
    checkValid(true);

    // CLEANUP -- Finally, fix up the surface mesh after octree operations
    SLIC_INFO("\tRegenerating the mesh");
    regenerateSurfaceMesh();

    SLIC_INFO("\tFinished generating the InOutOctree.");
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
bool InOutOctree<DIM>::hasGrayChild(const BlockIndex& blk) const
{
    for(int i=0; i < blk.numChildren(); ++i)
    {
        BlockIndex childBlk = blk.child(i);
        const InOutLeafData& childData = (*this)[childBlk];

        if(childData.hasTriangles())
            return true;
        else if(!childData.isLeaf() && hasGrayChild(childBlk) )
            return true;
    }

    return false;
}


template<int DIM>
bool InOutOctree<DIM>::colorLeafAndNeighbors(const BlockIndex& leafBlk, InOutLeafData& leafData)
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
                const InOutLeafData& neighborData = (*this)[neighborBlk];
                switch(neighborData.color())
                {
                case InOutLeafData::Black: leafData.setColor( InOutLeafData::Black ); break;
                case InOutLeafData::White: leafData.setColor( InOutLeafData::White ); break;
                case InOutLeafData::Gray:
                    if( withinGrayBlock(this->blockBoundingBox(leafBlk).centroid(),  neighborData) )
                        leafData.setColor(InOutLeafData::Black);
                    else
                        leafData.setColor(InOutLeafData::White);
                    break;
                case InOutLeafData::Undetermined:
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
                InOutLeafData& neighborData = (*this)[neighborBlk];
                if (! neighborData.isColored() )
                {
                    switch(leafData.color())
                    {
                    case InOutLeafData::Black:
                        neighborData.setColor(InOutLeafData::Black);
                        break;
                    case InOutLeafData::White:
                        neighborData.setColor(InOutLeafData::White);
                        break;
                    case InOutLeafData::Gray:
                        // PROBLEM HERE -- some neighbors are getting the wrong color!
                        if( withinGrayBlock(this->blockBoundingBox(neighborBlk).centroid(),  leafData) )
                            neighborData.setColor(InOutLeafData::Black);
                        else
                            neighborData.setColor(InOutLeafData::White);
                        break;
                    case InOutLeafData::Undetermined:
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
//
//            bool isGray = leafData.hasTriangles();
//
//            BlockIndex parent = BlockIndex(it->first, lev).parent();
//            bool hasGraySibling = hasGrayChild(parent);
//
//            SLIC_ASSERT( isGray || hasGraySibling );

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
bool InOutOctree<DIM>::withinGrayBlock(const SpacePt & pt, const InOutLeafData&  leafData) const
{
    SLIC_ASSERT( leafData.color() == InOutLeafData::Gray );
    SLIC_ASSERT( leafData.hasVertex() );

    const SpaceVector vec(vertexPosition( leafData.vertexIndex()), pt);

    // Iterate through triangles, count signs of dot products of vec with normals
    // If all dot products are positive, point is outside
    // If all dot products are negative, point is inside
    // Otherwise, short-circuit the calculation
    int sgnCount = 0;
    const int numTris = leafData.numTriangles();
    for(int i=0; i< numTris; ++i)
    {
        TriangleIndex tIdx = leafData.triangles()[i];

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

    TriangleIndex firstTriangle = leafData.triangles()[0];

    // Find a vertex from the local surface that is not incident in the first triangle
    TriVertIndices firstTriVerts = triangleVertexIndices( firstTriangle );
    TriVertIndices lastTriVerts = triangleVertexIndices( leafData.triangles()[ numTris -1]);
    VertexIndex otherVert = lastTriVerts[0];
    if( incidentInVertex(firstTriVerts, otherVert))
        otherVert = lastTriVerts[1];
    if( incidentInVertex(firstTriVerts, otherVert))
        otherVert = lastTriVerts[2];
    SLIC_ASSERT( !incidentInVertex(firstTriVerts, otherVert) );

    return ON_NEGATIVE_SIDE != orientation( vertexPosition(otherVert), trianglePositions( firstTriangle));
}


template<int DIM>
void InOutOctree<DIM>::insertMeshTriangles ()
{
    typedef asctoolkit::utilities::Timer Timer;

    typedef typename OctreeBaseType::MapType LeavesLevelMap;
    typedef typename LeavesLevelMap::iterator LeavesIterator;

    SLIC_ASSERT( !m_vertexSet.empty() );

    // Add all triangle references to the root
    BlockIndex rootBlock = this->root();
    InOutLeafData& rootData = (*this)[rootBlock];
    rootData.triangles().reserve( m_elementSet.size() );
    for(int idx=0; idx < m_elementSet.size(); ++idx)
        rootData.addTriangle(idx);

    // Iterate through octree levels and insert triangles into the blocks that they intersect
    for(int lev=0; lev< this->m_levels.size(); ++lev)
    {
        Timer levelTimer(true);

        LeavesLevelMap& levelLeafMap = this->m_leavesLevelMap[ lev ];
        for(LeavesIterator it=levelLeafMap.begin(), itEnd = levelLeafMap.end(); it != itEnd; ++it)
        {
            InOutLeafData& blkData = it->second;
            if(! blkData.hasTriangles() )
                continue;

            BlockIndex blk(it->first, lev);
            bool isInternal = !blkData.isLeaf();
            bool isLeafThatMustRefine = !isInternal && ! allTrianglesIncidentInCommonVertex(blk, blkData);

            // There is only work to do on internal blocks and leaves that must refine
            if( !isInternal && !isLeafThatMustRefine )
                continue;

            // At this point, we must distribute the triangles among the appropriate children
            InOutLeafData::TriangleList parentTriangles;
            parentTriangles.swap(blkData.triangles());
            int numTriangles = parentTriangles.size();

            // Refine the leaf if necessary
            if( isLeafThatMustRefine )
            {
                VertexIndex vIdx = blkData.vertexIndex();

                // Clear the current leaf and refine
                blkData.clear();
                this->refineLeaf(blk);

                // Note: At this stage, leaf might have triangles but no assigned vertex,
                //       so we need to check for that
                if( vIdx != InOutLeafData::NO_VERTEX && m_vertexToBlockMap[vIdx] == blk )
                    insertVertex(vIdx, blk.childLevel());
            }


            // Cache information about children blocks
            BlockIndex           childBlk[BlockIndex::NUM_CHILDREN];
            GeometricBoundingBox childBB[BlockIndex::NUM_CHILDREN];
            InOutLeafData*       childData[BlockIndex::NUM_CHILDREN];
            for(int j=0; j< BlockIndex::NUM_CHILDREN; ++j)
            {
                childBlk[j]  = blk.child(j);
                childBB[j]   = this->blockBoundingBox( childBlk[j] );
                childData[j] = &(*this)[childBlk[j] ];
				//childData[j]->reserveTriangles(numTriangles);
            }

            // Add all triangles to intersecting children blocks
            for(int i=0; i< numTriangles; ++i)
            {
                TriangleIndex tIdx = parentTriangles[i];
                SpaceTriangle spaceTri = trianglePositions(tIdx);
                GeometricBoundingBox tBB = triangleBoundingBox(tIdx);
                //tBB.scale(1.00001);

                for(int j=0; j< BlockIndex::numChildren(); ++j)
                {
                    if( childData[j]->isLeaf())
                    {
                        if(blockIndexesVertex(tIdx, childBlk[j]) || intersect(spaceTri, childBB[j]))
                            childData[j]->addTriangle(tIdx);
                    }
                    else if( intersect(tBB, childBB[j]) )
                    {
                        childData[j]->addTriangle(tIdx);
                    }

                }
            }

            blkData.clear();
        }

        if(! levelLeafMap.empty() )
          SLIC_DEBUG("\tInserting triangles into level " << lev <<" took " << levelTimer.elapsed() <<" seconds.");
    }

}


template<int DIM>
void InOutOctree<DIM>::insertVertex (VertexIndex idx, int startingLevel)
{
    static const double EPS_SQ = 1e-18;


    const SpacePt pt = getMeshVertexPosition(idx);

    BlockIndex block = this->findLeafBlock(pt, startingLevel);
    InOutLeafData& leafData = (*this)[block];

    #ifdef ATK_DEBUG
    if( idx == DEBUG_VERT_IDX)
    {
        SLIC_DEBUG("\t -- inserting pt with index " << idx << " and coords: " << pt
                  << ". Looking at block " << block
                  << " w/ blockBB " << this->blockBoundingBox(block)
                  << " indexing leaf vertex " << leafData.vertexIndex()
                    );
    }
    #endif

    if(! leafData.hasVertex())
    {
        leafData.setVertex(idx);

        // Update the vertex-to-block map for this vertex
        if(!m_vertexSet.empty())
            m_vertexToBlockMap[idx] = block;

    }
    else
    {
        // check if we should merge the vertices
        VertexIndex origVertInd = leafData.vertexIndex();
        if( squared_distance( pt, getMeshVertexPosition(origVertInd) ) >= EPS_SQ )
        {
            leafData.clear();
            this->refineLeaf(block);

            insertVertex(origVertInd, block.childLevel() );
            insertVertex(idx, block.childLevel() );
        }

    }

    #ifdef ATK_DEBUG
    if( leafData.vertexIndex() == DEBUG_VERT_IDX)
    {
        SLIC_DEBUG("-- vertex " << idx
                  << " is indexed in block " << block
                  << ". Leaf vertex is " << leafData.vertexIndex() );
    }
    #endif

}


template<int DIM>
void InOutOctree<DIM>::updateSurfaceMeshVertices()
{
    typedef asctoolkit::slam::Map<VertexIndex> IndexMap;

    // Create a map from old vertex indices to new vertex indices
    MeshVertexSet oldVerts( m_surfaceMesh->getMeshNumberOfNodes() );
    IndexMap vertexIndexMap( &oldVerts, InOutLeafData::NO_VERTEX);

    // Generate unique indices for new mesh vertices
    int uniqueVertexCounter = 0;
    for(int i=0; i< oldVerts.size(); ++i)
    {
        // Find the block and its indexed vertex in the octree
        BlockIndex leafBlock = this->findLeafBlock( getMeshVertexPosition(i) );
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
        SLIC_ASSERT( this->isLeaf(leafBlock) && (*this)[leafBlock].hasVertex() );

        (*this)[ leafBlock ].setVertex(i);
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
                                                         , InOutLeafData& leafData) const
{
    bool shareCommonVert = false;
    VertexIndex commonVert = InOutLeafData::NO_VERTEX;

    const int numTris = leafData.numTriangles();
    const InOutLeafData::TriangleList& tris = leafData.triangles();

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
            int sharedVert = InOutLeafData::NO_VERTEX;
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
        const InOutLeafData& data = (*this)[block];

        switch( data.color() )
        {
        case InOutLeafData::Black:      return true;
        case InOutLeafData::White:      return false;
        case InOutLeafData::Gray:       return withinGrayBlock( pt, data);
        case InOutLeafData::Undetermined:
            SLIC_ASSERT_MSG(false, "Error -- All leaf blocks must have a color. "
                        << " The color of leafBlock " << block
                        << " was 'Undetermined' when querying point " << pt  );
            break;
        }
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
            if(it->second.isLeaf())
                ++levelLeaves[lev];

            if(it->second.hasVertex())
                ++levelLeavesWithVert[ lev ];

            if(it->second.hasTriangles())
                levelTriangleRefCount[ lev ] += it->second.numTriangles();

            switch(it->second.color())
            {
            case InOutLeafData::Black:      ++levelBlackBlockCount[lev]; break;
            case InOutLeafData::White:      ++levelWhiteBlockCount[lev]; break;
            case InOutLeafData::Gray:       ++levelGrayBlockCount[lev];  break;
            case InOutLeafData::Undetermined:                            break;
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
    octreeStatsStr << "*** " << (trianglesAlreadyInserted ? "PM" : "PR") << " octree summary *** \n";

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
            if(trianglesAlreadyInserted)
            {
                octreeStatsStr << " Colors B,W,G ==> " << levelBlackBlockCount[lev]
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
                 << " (" <<   percentWithVert  << "% w/ vert); "
                 <<" \n\t There were " << totalTriangleRefCount << " triangle references "
                 <<" (avg. " << totalTriangleRefCount /  meshNumTriangles << " refs per triangle)."
                 <<"\n Colors B,W,G ==> " << totalBlackBlocks
                 << "," << totalWhiteBlocks
                 << "," << totalGrayBlocks
                 <<"\n";

    SLIC_INFO( octreeStatsStr.str() );

    if(trianglesAlreadyInserted)
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
                const InOutLeafData& leafData = it->second;
                if(leafData.isLeaf() && leafData.hasTriangles())
                {
                    const InOutLeafData::TriangleList& tris = leafData.triangles();
                    for(int i = 0; i < leafData.numTriangles(); ++i)
                    {
                        ++triCount[ tris[i]];

                        if(DEBUG_TRI_IDX == tris[i] )
                        {
                            BlockIndex blk(it->first, lev);
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
void InOutOctree<DIM>::checkValid(bool ATK_DEBUG_PARAM(trianglesAlreadyInserted)) const
{
#ifdef ATK_DEBUG
    typedef typename OctreeBaseType::MapType LeavesLevelMap;
    typedef typename LeavesLevelMap::const_iterator LeavesIterator;

    SLIC_DEBUG("Inside InOutOctree::checkValid() to verify state of PR or PM octree.");


    // Iterate through the tree
    // Internal blocks should not have associated vertex data
    // Leaf block consistency depends on 'color'
    //      Black or White blocks should have no vertex data (NO_VERTEX) and no triangles
    //      Gray blocks should have a vertex reference; it may or may not be located within the block
    //          The sum of vertices located within a block should equal the number of mesh vertices.
    //      Gray blocks should have one or more triangles.
    //          All triangles should be incident in a common vertex -- which equals the indexed vertex, if it exists.

    SLIC_ASSERT_MSG(!m_vertexSet.empty(), "InOutOctree::checkValid assumes that we have already inserted "
                    << "the vertices into the octree.");


    SLIC_DEBUG("--Checking that each vertex is in a leaf block of the tree.");
    const int numVertices = m_vertexSet.size();
    for(int i=0; i< numVertices; ++i)
    {
        // Check that we can find the leaf block indexing a vertex from its position
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

    // Check that each triangle is referenced in at least the three blocks indexing its vertices
    if(trianglesAlreadyInserted)
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
        SLIC_DEBUG("--Checking that internal blocks have no data, and that leaves satisfy all PM conditions");
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

    SLIC_DEBUG("done.");

#endif
}



} // end namespace quest

#endif  // SPATIAL_OCTREE__HXX_
