/*
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

#ifndef QUEST_IMPLICIT_GRID__HPP_
#define QUEST_IMPLICIT_GRID__HPP_

#include "axom/config.hpp"
#include "axom_utils/Utilities.hpp"  // for clamp functions

#include "slic/slic.hpp"

#include "primal/BoundingBox.hpp"
#include "primal/Point.hpp"
#include "primal/Vector.hpp"
#include "primal/RectangularLattice.hpp"

// #include "slam/SizePolicies.hpp"
#include "slam/OrderedSet.hpp"
#include "slam/Map.hpp"

#include <vector>

#ifdef AXOM_USE_BOOST
#  include <boost/dynamic_bitset.hpp>
#else
#  error quest::ImplicitGrid uses boost headers.
#endif

namespace axom {
namespace quest {


/*!
 * \class ImplicitGrid
 *
 * \brief An implicit grid is an occupancy-based spatial index over an indexed set
 * of objects in space.
 *
 * An ImplicitGrid divides a given rectilinear slab of space (defined by an
 * axis aligned bounding box) into uniformly sized cells of a specified resolution.
 * The GridCells of the ImplicitGrid index a subset of the items from an indexed
 * set (whose cardinality is specified during the ImplicitGrid's initialization).
 * Users can insert items from the indexed set into an ImplicitGrid by providing
 * the item's bounding box and index.
 *
 * In contrast to a primal::UniformGrid, which encodes an array of indices
 * for each cell in the underlying multidimensional grid,
 * the ImplicitGrid encodes a single array of bins per dimension, each of which
 * has a bitset over the index space.  Thus, the storage overhead of an
 * ImplicitGrid is fixed at initialization time to numElts * sum_i ( res[i] ) bits.
 * Queries are implemented in terms of unions and intersections of bitsets.
 *
 * One might prefer an ImplicitGrid over a UniformGrid when one expects
 * a relatively dense index relative to the grid resolution (i.e. that
 * there will be many items indexed per bucket).  The ImplicitGrid
 * is designed for quick indexing and searching over a static (and relatively
 * small index space) in a relatively coarse grid.
 */
template<int NDIMS, typename TheIndexType = int>
class ImplicitGrid
{
public:
  typedef TheIndexType                          IndexType;
  typedef primal::Point<IndexType, NDIMS>       GridCell;
  typedef primal::Point<double, NDIMS>          SpacePoint;
  typedef primal::Vector<double, NDIMS>         SpaceVec;

  typedef primal::BoundingBox<double, NDIMS>    SpatialBoundingBox;
  typedef primal::RectangularLattice<NDIMS, double, IndexType> LatticeType;

  typedef slam::OrderedSet< slam::policies::RuntimeSize<IndexType> > ElementSet;
  typedef slam::OrderedSet< slam::policies::RuntimeSize<IndexType> > BinSet;


  typedef boost::dynamic_bitset<>  BitsetType;
  typedef slam::Map<BitsetType>    BinBitMap;

  /*!
   * \brief Default constructor for an ImplicitGrid
   *
   * \note Users must call initialize() to initialize the ImplicitGrid
   *       after constructing with the default constructor
   */
  ImplicitGrid(): m_initialized(false) {}

  /*!
   * \brief Constructor for an implicit grid from a bounding box, a grid resolution a number of elements
   */
  ImplicitGrid(const SpatialBoundingBox& boundingBox, const GridCell& gridRes, int numElts)
    : m_bb(boundingBox),
      m_lattice( primal::rectangular_lattice_from_bounding_box(boundingBox, gridRes.array()))
  {
     m_elementSet = ElementSet(numElts);
     for(int i=0; i<NDIMS; ++i)
     {
        m_bins[i] = BinSet(gridRes[i]);
        m_binData[i] = BinBitMap(&m_bins[i], BitsetType(m_elementSet.size()));
     }

    // Set the expansion factor for each element to a small fraction of the grid's bounding boxes diameter
    // TODO: Add a constructor that allows users to set the expansion factor
    static const double EPS = 1e-8;
    m_expansionFactor = m_bb.range().norm() * EPS;

    m_initialized = true;
  }
  
  /*! Predicate to check if the ImplicitGrid has been initialized */
  bool isInitialized() const { return m_initialized; }


  /*!
   * Initializes an implicit grid or resolution gridRes over an axis aligned domain covered by boundingBox.
   * The implicit grid indexes a set with numElts elements.
   *
   * \pre The ImplicitGrid has not already been initialized
   */
  void initialize(const SpatialBoundingBox& boundingBox, const GridCell& gridRes, int numElts) 
  {
    SLIC_ASSERT( !m_initialized);

    m_bb = boundingBox;
    m_lattice = primal::rectangular_lattice_from_bounding_box(boundingBox, gridRes.array());
    m_elementSet = ElementSet(numElts);

    for(int i=0; i<NDIMS; ++i)
    {
      m_bins[i] = BinSet(gridRes[i]);
      m_binData[i] = BinBitMap(&m_bins[i], BitsetType(m_elementSet.size()));
    }

    // Set the expansion factor for each element to a small fraction of the grid's bounding boxes diameter
    // TODO: Add a constructor that allows users to set the expansion factor
    static const double EPS = 1e-8;
    m_expansionFactor = m_bb.range().norm() * EPS;

    m_initialized = true;
  }

  /*!
   * \brief Inserts an element with index id and bounding box bbox into the implicit grid
   *
   * \note bbox is intentionally passed by value since insert() modifies its bounds
   */
  void insert(SpatialBoundingBox bbox, int id)
  {
    SLIC_ASSERT(m_initialized);
    
    // Note: We slightly inflate the bbox of the objects.
    //       This effectively ensures that objects on grid boundaries are added all nearby grid cells.

    bbox.expand(m_expansionFactor);
    // SLIC_INFO("Inserting element " << id << " with bb " << bbox);


    const GridCell lowerCell = m_lattice.gridCell( bbox.getMin() );
    const GridCell upperCell = m_lattice.gridCell( bbox.getMax() );

    for(int i=0; i< NDIMS; ++i)
    {
      BinBitMap& binData = m_binData[i];

      const IndexType lower = axom::utilities::clampLower(lowerCell[i], IndexType() );
      const IndexType upper = axom::utilities::clampUpper(upperCell[i], highestBin(i) );
      for(int j= lower; j <= upper; ++j)
      {
        binData[j].set(id);
      }
    }
  }

  /*! Finds the candidate elements in the vicinity of query point pt  */
  BitsetType getCandidates(const SpacePoint& pt) const
  {
    if(! m_bb.contains(pt) )
      return BitsetType(0);

    const GridCell gridCell = m_lattice.gridCell(pt);

    // Note: Need to clamp the upper range of the gridCell
    //       to handle points on the upper boundaries of the bbox
    //       This is valid since we've already ensured that pt is in the bbox.
    IndexType idx = axom::utilities::clampUpper(gridCell[0], highestBin(0));
    BitsetType res = m_binData[0][ idx ];

    for(int i=1; i< NDIMS; ++i)
    {
      idx = axom::utilities::clampUpper(gridCell[i], highestBin(i));
      res &= m_binData[i][idx];
    }

    return res;
  }

  /*! Returns the list of candidates as an explicit list of IndexType */
  std::vector<IndexType> getCandidatesAsArray(const SpacePoint& pt) const
  {
    std::vector<IndexType> candidatesVec;

    BitsetType candidateBits = getCandidates(pt);
    candidatesVec.reserve( candidateBits.count() );
    for(std::size_t eltIdx = candidateBits.find_first(); eltIdx != BitsetType::npos; eltIdx = candidateBits.find_next( eltIdx) )
    {
      candidatesVec.push_back( eltIdx );
    }

    return candidatesVec;
  }


  /*!
   * Tests whether grid cell gridPt indexes the element with index id
   */
  bool contains(const GridCell& gridCell, int id) const
  {
    bool ret = true;

    if(! m_elementSet.isValidIndex(id) )
      ret = false;

    for(int i=0; i< NDIMS; ++i)
    {
      ret = ret
          && m_bins[i].isValidIndex(gridCell[i])
          && m_binData[ i][ gridCell[i] ].test( id);
    }

    return ret;
  }

private:

  /*!
   * \brief Returns the bin index in the given dimension dim
   *
   * \pre 0 <= dim < NDIMS
   */
  IndexType highestBin(int dim) const
  {
    SLIC_ASSERT(0 <= dim && dim < NDIMS);
    return m_bins[dim].size()-1;
  }

private:

  SpatialBoundingBox m_bb; //!< The bounding box of the ImplicitGrid
  LatticeType        m_lattice; //!< A lattice to help in converting from points in space to GridCells
  double             m_expansionFactor; //!< The amount by which to expand bounding boxes

  ElementSet         m_elementSet; //!< The index set of the elements
  BinSet             m_bins[NDIMS]; //!< A set of bins, per dimension
  BinBitMap          m_binData[NDIMS]; //!< The data associated with each bin

  bool               m_initialized; //!< Tracks whether the ImplicitGrid has been initialized

};



} // end namespace quest
} // end namespace axom

#endif  // QUEST_IMPLICIT_GRID__HPP_
