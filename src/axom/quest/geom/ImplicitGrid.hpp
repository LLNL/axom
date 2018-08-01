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

#ifndef QUEST_IMPLICIT_GRID__HPP_
#define QUEST_IMPLICIT_GRID__HPP_

#include "axom/config.hpp"
#include "axom/core/Types.hpp"            // AXOM_NULLPTR
#include "axom/core/utilities/Utilities.hpp"  // for clamp functions

#include "axom/slic/interface/slic.hpp"

#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/RectangularLattice.hpp"

#include "axom/slam/policies/SizePolicies.hpp"
#include "axom/slam/OrderedSet.hpp"
#include "axom/slam/Map.hpp"
#include "axom/slam/BitSet.hpp"

#include <vector>

namespace axom
{
namespace quest
{


/*!
 * \class ImplicitGrid
 *
 * \brief An implicit grid is an occupancy-based spatial index over an indexed
 * set of objects in space.
 *
 * An ImplicitGrid divides a given rectilinear slab of space (defined by an
 * axis aligned bounding box) into uniformly sized cells
 * of a specified resolution.
 * The GridCells of the ImplicitGrid index a subset of the items from an indexed
 * set whose cardinality is specified during the ImplicitGrid's initialization.
 * Users can insert items from the indexed set into an ImplicitGrid by providing
 * the item's bounding box and index.
 *
 * In contrast to a primal::UniformGrid, which encodes an array of indices
 * for each cell in the underlying multidimensional grid,
 * the ImplicitGrid encodes a single array of bins per dimension, each of which
 * has a bitset over the index space.  Thus, the storage overhead of an
 * ImplicitGrid is fixed at initialization time to
 *   \f$ numElts * sum_i { res[i] } \f$ bits.
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
  typedef TheIndexType IndexType;
  typedef primal::Point<IndexType, NDIMS>       GridCell;
  typedef primal::Point<double, NDIMS>          SpacePoint;
  typedef primal::Vector<double, NDIMS>         SpaceVec;

  typedef primal::BoundingBox<double, NDIMS>    SpatialBoundingBox;
  typedef primal::RectangularLattice<NDIMS, double, IndexType> LatticeType;

  typedef slam::policies::RuntimeSize<IndexType> SizePolicy;
  typedef slam::OrderedSet< SizePolicy > ElementSet;
  typedef slam::OrderedSet< SizePolicy > BinSet;


  typedef slam::BitSet BitsetType;
  typedef slam::Map<BitsetType>    BinBitMap;

  /*!
   * \brief Default constructor for an ImplicitGrid
   *
   * \note Users must call initialize() to initialize the ImplicitGrid
   *       after constructing with the default constructor
   */
  ImplicitGrid() : m_initialized(false) {}

  /*!
   * \brief Constructor for an implicit grid from a bounding box,
   * a grid resolution a number of elements.
   *
   * \param [in] boundingBox Bounding box of domain to index
   * \param [in] gridRes Pointer to resolution for lattice covering bounding box
   * \param [in] numElts The number of elements to be indexed
   *
   * \pre \a gridRes is either NULL or has \a NDIMS coordinates
   * \sa initialize() for details on setting grid resolution
   * when \a gridRes is NULL
   */
  ImplicitGrid(const SpatialBoundingBox& boundingBox,
               const GridCell* gridRes,
               int numElts)
    : m_bb(boundingBox), m_initialized(false)
  {
    initialize(m_bb, gridRes, numElts);
  }

  /*!
   * \brief Constructor for an implicit grid from arrays of primitive types
   *
   * \param [in] bbMin Lower bounds of mesh bounding box
   * \param [in] bbMax Upper bounds of mesh bounding box
   * \param [in] gridRes Resolution for lattice covering mesh bounding box
   * \param [in] numElts The number of elements in the index space
   *
   * \pre \a bbMin and \a bbMax are not NULL and have \a NDIMS coordinates
   * \pre \a gridRes is either NULL or has \a NDIMS coordinates
   * \sa initialize() for details on setting grid resolution
   * when \a gridRes is NULL
   */
  ImplicitGrid(const double* bbMin,
               const double* bbMax,
               const int* gridRes,
               int numElts)
    : m_initialized(false)
  {
    SLIC_ASSERT( bbMin != nullptr);
    SLIC_ASSERT( bbMax != nullptr);

    // Set up the grid resolution from the gridRes array
    //   if NULL, GridCell parameter to initialize should also be NULL
    const GridCell* pRes =
      (gridRes != nullptr) ? &m_gridRes : nullptr;

    initialize(
      SpatialBoundingBox( SpacePoint(bbMin), SpacePoint(bbMax) ),
      pRes, numElts);
  }

  /*! Predicate to check if the ImplicitGrid has been initialized */
  bool isInitialized() const { return m_initialized; }


  /*!
   * \brief Initializes an implicit grid or resolution gridRes over an axis
   * aligned domain covered by boundingBox. The implicit grid indexes a set
   * with numElts elements.
   *
   * \param [in] boundingBox Bounding box of domain to index
   * \param [in] gridRes Resolution for lattice covering bounding box
   * \param [in] numElts The number of elements to be indexed
   * \pre The ImplicitGrid has not already been initialized
   *
   * \note When \a gridRes is NULL, we use a heuristic to set the grid
   * resolution to the N^th root of \a numElts. We also ensure that
   * the resolution along any dimension is at least one.
   *
   */
  void initialize(const SpatialBoundingBox& boundingBox,
                  const GridCell* gridRes,
                  int numElts)
  {
    SLIC_ASSERT( !m_initialized);

    // Setup Grid Resolution, dealing with possible null pointer
    if(gridRes == nullptr)
    {
      // Heuristic: use the n^th root of the number of elements
      // add 0.5 to round to nearest integer
      double nthRoot = std::pow(static_cast<double>(numElts), 1./ NDIMS );
      int dimRes = static_cast<int>(nthRoot + 0.5);
      m_gridRes = GridCell(dimRes);
    }
    else
    {
      m_gridRes = GridCell(*gridRes);
    }

    // ensure that resolution in each dimension is at least one
    for(int i=0 ; i< NDIMS ; ++i)
    {
      m_gridRes[i] = axom::utilities::max( m_gridRes[i], 1);
    }

    // Setup lattice
    m_bb = boundingBox;
    m_lattice = primal::rectangular_lattice_from_bounding_box(boundingBox,
                                                              m_gridRes.array());
    m_elementSet = ElementSet(numElts);

    for(int i=0 ; i<NDIMS ; ++i)
    {
      m_bins[i] = BinSet(m_gridRes[i]);
      m_binData[i] = BinBitMap(&m_bins[i], BitsetType(m_elementSet.size()));
    }

    // Set the expansion factor for each element to a small fraction of the
    // grid's bounding boxes diameter
    // TODO: Add a constructor that allows users to set the expansion factor
    const double EPS = 1e-8;
    m_expansionFactor = m_bb.range().norm() * EPS;

    m_initialized = true;
  }

  /*! Accessor for ImplicitGrid's resolution */
  const GridCell& gridResolution() const { return m_gridRes; }

  /*! Returns the number of elements in the ImplicitGrid's index set */
  int numIndexElements() const { return m_elementSet.size(); }

  /*!
   * \brief Inserts an element with index \a idx and bounding box \a bbox
   * into the implicit grid
   *
   * \param [in] bbox The bounding box of the element
   * \param [in] idx  The index of the element
   *
   * \note \a bbox is intentionally passed by value since insert()
   * modifies its bounds
   */
  void insert(SpatialBoundingBox bbox, IndexType idx)
  {
    SLIC_ASSERT(m_initialized);

    // Note: We slightly inflate the bbox of the objects.
    //       This effectively ensures that objects on grid boundaries are added
    //       in all nearby grid cells.

    bbox.expand(m_expansionFactor);

    const GridCell lowerCell = m_lattice.gridCell( bbox.getMin() );
    const GridCell upperCell = m_lattice.gridCell( bbox.getMax() );

    for(int i=0 ; i< NDIMS ; ++i)
    {
      BinBitMap& binData = m_binData[i];

      const IndexType lower =
        axom::utilities::clampLower(lowerCell[i], IndexType() );
      const IndexType upper =
        axom::utilities::clampUpper(upperCell[i], highestBin(i) );

      for(int j= lower ; j <= upper ; ++j)
      {
        binData[j].set(idx);
      }
    }
  }

  /*!
   * Finds the candidate elements in the vicinity of query point pt
   *
   * \param [in] pt The query point
   * \return A bitset \a bSet whose bits correspond to
   * the elements of the IndexSet.
   * The bits of \a bSet are set if their corresponding element bounding boxes
   * overlap the grid cell containing \a pt.
   */
  BitsetType getCandidates(const SpacePoint& pt) const
  {
    if(!m_initialized || !m_bb.contains(pt) )
      return BitsetType(0);

    const GridCell gridCell = m_lattice.gridCell(pt);

    // Note: Need to clamp the upper range of the gridCell
    //       to handle points on the upper boundaries of the bbox
    //       This is valid since we've already ensured that pt is in the bbox.
    IndexType idx = axom::utilities::clampUpper(gridCell[0], highestBin(0));
    BitsetType res = m_binData[0][ idx ];

    for(int i=1 ; i< NDIMS ; ++i)
    {
      idx = axom::utilities::clampUpper(gridCell[i], highestBin(i));
      res &= m_binData[i][idx];
    }

    return res;
  }

  /*!
   * Returns the list of candidates as an explicit list of IndexType
   *
   * \param [in] pt The query point
   * \return An list of indexes from the IndexSet whose corresponding
   * bounding boxes overlap the grid cell containing \a pt.
   *
   * \note This function returns the same indices as \a getCandidates()
   * But the results here are converted into an explicit list.
   */
  std::vector<IndexType> getCandidatesAsArray(const SpacePoint& pt) const
  {
    std::vector<IndexType> candidatesVec;

    BitsetType candidateBits = getCandidates(pt);
    candidatesVec.reserve( candidateBits.count() );
    for(IndexType eltIdx = candidateBits.find_first() ;
        eltIdx != BitsetType::npos ;
        eltIdx = candidateBits.find_next( eltIdx) )
    {
      candidatesVec.push_back( eltIdx );
    }

    return candidatesVec;
  }


  /*!
   * Tests whether grid cell gridPt indexes the element with index idx
   *
   * \param [in] gridCell The cell within the ImplicitGrid that we are testing
   * \param [in] idx An element index from the IndexSet to test
   *
   * \pre \a idx should be in the IndexSet.  I.e. 0 <= idx < numIndexElements()
   * \return True if the bounding box of element \a idx overlaps
   * with GridCell \a gridCell.
   */
  bool contains(const GridCell& gridCell, IndexType idx) const
  {
    bool ret = true;

    if(!m_elementSet.isValidIndex(idx) )
      ret = false;

    for(int i=0 ; i< NDIMS ; ++i)
    {
      ret = ret
            && m_bins[i].isValidIndex(gridCell[i])
            && m_binData[ i][ gridCell[i] ].test( idx);
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

  //! The bounding box of the ImplicitGrid
  SpatialBoundingBox m_bb;

  //! A lattice to help in converting from points in space to GridCells
  LatticeType m_lattice;

  //! The amount by which to expand bounding boxes
  double m_expansionFactor;

  //! Resolution of the ImplicitGrid
  GridCell m_gridRes;

  //! The index set of the elements
  ElementSet m_elementSet;

  //! A set of bins, per dimension
  BinSet m_bins[NDIMS];

  //! The data associated with each bin
  BinBitMap m_binData[NDIMS];

  //! Tracks whether the ImplicitGrid has been initialized
  bool m_initialized;
};


} // end namespace quest
} // end namespace axom

#endif  // QUEST_IMPLICIT_GRID__HPP_
