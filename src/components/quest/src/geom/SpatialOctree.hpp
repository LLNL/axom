#ifndef SPATIAL_OCTREE__HXX_
#define SPATIAL_OCTREE__HXX_

#include "quest/BoundingBox.hpp"
#include "quest/Point.hpp"
#include "quest/Vector.hpp"

#include "slic/slic.hpp"

#include "slam/Map.hpp"

#include "quest/Mesh.hpp"
#include "quest/OctreeBase.hpp"


namespace quest
{



/**
 * \class
 * \brief Adds spatial extents to an OctreeBase, allowing point location
 */
template<int DIM, typename LeafNodeType>
class SpatialOctree : public OctreeBase<DIM, LeafNodeType>
{
public:

    typedef quest::BoundingBox<double,DIM> GeometricBoundingBox;
    typedef quest::Point<double,DIM> SpacePt;
    typedef quest::Vector<double,DIM> SpaceVector;

    typedef typename OctreeBase<DIM, LeafNodeType>::GridPt GridPt;
    typedef typename OctreeBase<DIM, LeafNodeType>::MapType MapType;
    typedef typename OctreeBase<DIM, LeafNodeType>::BlockIndex BlockIndex;

//    typedef typename OctreeBase<DIM, LeafNodeType>::OctreeLevels OctreeLevels;

    typedef asctoolkit::slam::Map<SpaceVector> SpaceVectorLevelMap;

public:
    /**
     * \brief Construct a spatial octree from a spatial bounding box
     * \param [in] bb The spatial extent to be indexed by the octree
     */
    SpatialOctree(const GeometricBoundingBox& bb)
        : OctreeBase<DIM,LeafNodeType>()
        , m_deltaLevelMap(& this->m_levels)
        , m_boundingBox(bb)
    {
        SpaceVector bbRange = bb.range();
        for(int lev = 0; lev < this->m_levels.size(); ++lev)
        {
            m_deltaLevelMap[lev] = bbRange / static_cast<double>(1<<lev);
        }

    }

    /**
     * \brief Return the spatial bounding box of a grid cell at the given level or resolution
     */
    GeometricBoundingBox blockBoundingBox(const BlockIndex & block)
    {
        return blockBoundingBox( block.pt(), block.level() );
    }

    /**
     * \brief Return the spatial bounding box of a grid cell at the given level or resolution
     */
    GeometricBoundingBox blockBoundingBox(const GridPt & gridPt, int level)
    {
        const SpaceVector& deltaVec = m_deltaLevelMap[ level];

        SpacePt lower(m_boundingBox.getMin());
        SpacePt upper(m_boundingBox.getMin());
        for(int i=0; i< DIM; ++i)
        {
            lower[i] +=  gridPt[i]   * deltaVec[i];
            upper[i] += (gridPt[i]+1)* deltaVec[i];
        }

        return GeometricBoundingBox( lower,upper );
    }


    /**
     * Returns the width of an octree block at level of resolution level
     */
    const SpaceVector& spacingAtLevel(int level) const
    {
        return m_deltaLevelMap[level];
    }

    /**
     * \brief Finds the index of the leaf block covering the query point pt
     * \param [in] pt The query point in space
     * \param [in] startingLevel Optional starting level for the query
     * \pre pt must be in the bounding box of the octree (i.e. boundingBox.contains(pt) == true )
     * \note The collection of leaves covers the bounding box, and the interiors of the leaves do not
     * intersect, so every point in the bounding box should be located in a unique leaf block.
     * \note We are assuming a half-open interval on the bounding boxes.
     * \return The block index (i.e. grid point and level) of the leaf block containing the query point
     */
    BlockIndex findLeafBlock(const SpacePt& pt, int startingLevel = 0) const
    {
        bool found = false;
        BlockIndex leafBlock = BlockIndex::invalid_index();

        SLIC_ASSERT( m_boundingBox.contains(pt) );

        for(int lev=startingLevel; !found && lev < this->m_levels.size(); ++lev)
        {
            GridPt gridPt = findGridCellAtLevel(pt, lev);
            found = this->isLeaf(gridPt, lev);

            if(found)
                leafBlock = BlockIndex(gridPt, lev);
        }

        SLIC_ASSERT_MSG(found, "Point " << pt << " not found "
                        <<"in a leaf block of the octree");

        return leafBlock;
    }

    /**
     * \brief Utility function to find the quantized level lev grid cell of Point pt
     * \param [in] pt The point at which we are querying.
     * \param [in] level The level or resolution.
     * \pre \f$ 0 \le lev < octree.maxLeafLevel() \f$
     * \return The grid point of the block covering this point at this level
     * \todo KW: Should this function be protected? Is it generally useful?
     */
    GridPt findGridCellAtLevel(const SpacePt& pt, int level) const
    {
        SpaceVector ptVec(m_boundingBox.getMin(), pt);
        const SpaceVector& deltaVec = m_deltaLevelMap[ level];

        // Find the octree block that covers us at this level
        // Note: we are assuming a half-open interval for all blocks
        //       that are not on the upper boundary of the domain
        GridPt quantizedPt = elementwiseQuantizedRatio( ptVec, deltaVec);

        // Fix for when we are on the upper boundary of the domain
        for(int i=0; i< DIM; ++i)
        {
            if( quantizedPt[i] == 1<<level )
                --quantizedPt[i];
        }

        return quantizedPt;
    }

private:
    /**
     * \brief Helper function to quantize to the integer grid
     */
    GridPt elementwiseQuantizedRatio(const SpaceVector& ptFromBBMin, const SpaceVector&  cellWidth) const
    {
        GridPt gridPt;
        for(int i=0; i< DIM; ++i)
            gridPt[i] = std::floor( ptFromBBMin[i] / cellWidth[i]);

        return gridPt;
    }

private:
  DISABLE_COPY_AND_ASSIGNMENT(SpatialOctree)

protected:
    SpaceVectorLevelMap     m_deltaLevelMap;
    GeometricBoundingBox    m_boundingBox;
};


} // end namespace quest

#endif  // SPATIAL_OCTREE__HXX_
