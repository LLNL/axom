#ifndef INOUT_OCTREE__HXX_
#define INOUT_OCTREE__HXX_

#include "quest/BoundingBox.hpp"
#include "quest/Point.hpp"
#include "quest/Vector.hpp"
#include "quest/SquaredDistance.hpp"
#include "quest/fuzzy_compare.hpp"


#include "slic/slic.hpp"

#include "slam/Map.hpp"

#include "quest/Mesh.hpp"
#include "quest/SpatialOctree.hpp"


#ifdef USE_CXX11
  // Note: Replace the explicit timer functionality once
  //       a timer class is added to common.
  #include <chrono>
  #include <ctime>
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
        enum { NO_VERTEX = -1 };
        typedef int VertexIndex;

    public:
        /**
         * \brief Constructor for an InOutOctree LeafData
         * \param vInd The index of a vertex (optional; default is to not set a vertex)
         */
        InOutLeafData(VertexIndex vInd = NO_VERTEX) : LeafData(), m_vertIndex(vInd) {}

        /**
         * \brief Copy constructor for an InOutLeafData instance
         */
        InOutLeafData(const InOutLeafData& other) : LeafData(other), m_vertIndex(other.m_vertIndex) {}

        /**
         * \brief Assignment operator for an InOutLeafData instance
         */
        InOutLeafData& operator=(const InOutLeafData& other)
        {
            LeafData::operator=(other);
            this->m_vertIndex = other.m_vertIndex;
            return *this;
        }

        /**
         * \brief Checks whether there is a vertex associated with this leaf
         */
        bool hasVertex() const { return m_vertIndex != NO_VERTEX; }

        /**
         * \brief Removes all indexed data from this leaf
         */
        void clear() { m_vertIndex = NO_VERTEX; }

        /**
         * \brief Sets the vertex associated with this leaf
         */
        void setVertex(VertexIndex vInd) { m_vertIndex = vInd; }

        /**
         * \brief Accessor for the index of the vertex associated with this leaf
         */
        VertexIndex& vertexIndex() { return m_vertIndex; }

        /**
         * \brief Constant accessor for the index of the vertex associated with this leaf
         */
        const VertexIndex& vertexIndex() const { return m_vertIndex; }

        /**
         * \brief Equality operator to determine if two InOutLeafData instances are equivalent
         */
        friend bool operator==(const InOutLeafData& lhs, const InOutLeafData& rhs )
        {
            return (static_cast<const LeafData&>(lhs) == static_cast<const LeafData&>(rhs))
                && (lhs.m_vertIndex == rhs.m_vertIndex);
        }
    private:

        VertexIndex m_vertIndex;
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


    typedef meshtk::Mesh SurfaceMesh;
    typedef int VertexIndex;
    typedef int TriangleIndex;


public:
    /**
     * \brief Construct an InOutOctree to handle containment queries on a surface mesh
     * \param [in] bb The spatial extent covered by the octree
     */
    InOutOctree(const GeometricBoundingBox& bb, SurfaceMesh* meshPtr)
        : SpatialOctreeType(bb)
        , m_surfaceMesh(meshPtr)
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
     */
    void insertVertex(VertexIndex idx);

    /**
     * \brief Helper function to insert a triangle into the octree
     * \param idx The index of the triangle that we are inserting
     */
    void insertTriangle(TriangleIndex idx);

    /**
     * \brief Helper function to retrieve the position of the vertex from the mesh
     */
    SpacePt vertexPosition(VertexIndex idx) const
    {
        SpacePt pt;
        m_surfaceMesh->getMeshNode(idx, pt.data() );
        return pt;
    }

    /**
     * Utility function to print some statistics about the InOutOctree instance
     */
    void printOctreeStats() const;

private:
  DISABLE_COPY_AND_ASSIGNMENT(InOutOctree)

protected:
    SurfaceMesh* m_surfaceMesh;

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

    int numMeshVerts = m_surfaceMesh->getMeshNumberOfNodes();
    for(int idx=0; idx < numMeshVerts; ++idx)
    {
        insertVertex(idx);
    }

#ifdef USE_CXX11
    vEnd = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = vEnd-vStart;
    SLIC_INFO("\tInserting vertices took " << elapsed_seconds.count() << " seconds.");
#endif

    printOctreeStats();
}


template<int DIM>
void InOutOctree<DIM>::insertVertex (VertexIndex idx)
{
    static const double EPS = 1e-18;


    SpacePt pt = vertexPosition(idx);

    BlockIndex block = this->findLeafBlock(pt);
    InOutLeafData& leafData = (*this)[block];

    if(! leafData.hasVertex())
    {
        leafData.setVertex(idx);
    }
    else
    {
        // check if we should merge the vertices
        VertexIndex origVertInd = leafData.vertexIndex();
        SpacePt otherPt = vertexPosition(origVertInd);

        if( squared_distance( pt, otherPt) >= EPS )
        {
            leafData.clear();

            this->refineLeaf(block);
            insertVertex(origVertInd);
            insertVertex(idx);
        }

    }

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

    int totalBlocks = 0;
    int totalLeaves = 0;
    int totalLeavesWithVert = 0;

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
        }

        totalBlocks += levelBlocks[lev];
        totalLeaves += levelLeaves[lev];
        totalLeavesWithVert += levelLeavesWithVert[lev];
    }


    std::stringstream octreeStatsStr;
    octreeStatsStr << "*** PR octree summary *** \n";
    for(int lev=0; lev< this->m_levels.size(); ++lev)
    {
        if(levelBlocks[lev] > 0)
        {
            octreeStatsStr << "\t Level " << lev
                           << " has " << levelLeavesWithVert[lev] << " leaves with vert"
                           << " out of " << levelLeaves[lev] << " leaves;"
                           << " and " << levelBlocks[lev] - levelLeaves[lev] << " internal blocks.\n";
        }
    }
    octreeStatsStr<<"  Mesh has " << m_surfaceMesh->getMeshNumberOfNodes() << " vertices."
                 <<"\n  Octree has " << totalLeavesWithVert << " filled leaves; "
                 <<  totalLeaves << " leaves; "
                 <<  totalBlocks - totalLeaves << " internal blocks"
                 <<" and " << totalBlocks << " overall blocks.\n";

    SLIC_INFO( octreeStatsStr.str() );
}


template<int DIM>
void InOutOctree<DIM>::insertTriangle (TriangleIndex idx)
{
    SLIC_ERROR("Not implemented yet.");
}

template<int DIM>
bool InOutOctree<DIM>::within(const SpacePt& pt) const
{
    if( ! blockBoundingBox( this->root() ).contains(pt) )
        return false;

    SLIC_ERROR("Not implemented yet.");
}


} // end namespace quest

#endif  // SPATIAL_OCTREE__HXX_
