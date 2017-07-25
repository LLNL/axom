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

#include "slic/slic.hpp"

#include "primal/BoundingBox.hpp"
#include "primal/Point.hpp"
#include "primal/Vector.hpp"
#include "primal/RectangularLattice.hpp"

// #include "slam/SizePolicies.hpp"
#include "slam/OrderedSet.hpp"
#include "slam/Map.hpp"

// #include <ostream>   // for ostream in print
#include <cmath>  // for std::floor
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
 * \brief An implicit grid is an occupancy-based spatial index over a set of indices
 */
template<int NDIMS, typename IndexType = int>
class ImplicitGrid
{
public:
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
   * \brief Default constructor for an implicit grid
   *
   * \post The ImplicitGrid is not initialized after the default constructor
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
  
  /*! Predicate to see if the implicit grid has been initialized */
  bool isInitialized() const { return m_initialized; }


  /*!
   * Initializes an implicit grid or resolution gridRes over an axis aligned domain covered by boundinbBox.
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
   * Inserts an element with index id and bounding box bbox into the implicit grid
   * \note bbox is intentionally passed by value since we will slightly inflate its bounds
   */
  void insert(SpatialBoundingBox bbox, int id)
  {
    SLIC_ASSERT(m_initialized);
    
    // Note: We slightly inflate the bbox of the objects.
    //       This effectively ensures that objects on grid boundaries are added all nearby grid cells.

    bbox.expand(m_expansionFactor);
    // SLIC_INFO("Inserting element " << id << " with bb " << bbox);


    GridCell lowerCell = m_lattice.gridCell( bbox.getMin() );
    GridCell upperCell = m_lattice.gridCell( bbox.getMax() );

    for(int i=0; i< NDIMS; ++i)
    {
      BinBitMap& binData = m_binData[i];
      IndexType lower = clampLower(lowerCell[i], m_bins[i] );
      IndexType upper = clampUpper(upperCell[i], m_bins[i] );

      //SLIC_INFO("Updated bounds l0 -- hi -> " << lower << " -- " << upper);
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

    GridCell gridCell = m_lattice.gridCell(pt);

    // Note: Need to clamp the upper range of the gridCell
    //       to handle points on the upper boundaries of the bbox
    //       This is valid since we've already ensured that pt is in the bbox.
    BitsetType res = m_binData[0][ clampUpper( gridCell[0], m_bins[0]) ];

    for(int i=1; i< NDIMS; ++i)
    {
       res &= m_binData[i][ clampUpper( gridCell[i], m_bins[i]) ];
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
   * \brief Utility function to ensure that a space point is within a valid bin
   *
   * \param val The current bin index
   * \param binSet The set of valid bins
   * \return The min of val and binSet.size()-1
   */
  IndexType clampUpper(IndexType val, const BinSet& binSet) const
  {
    const IndexType upperVal = binSet.size() -1;
    return val > upperVal ? upperVal : val;
  }

  /*!
   * \brief Utility function to ensure that a space point is within a valid bin
   *
   * \param val The current bin index
   * \param binSet The set of valid bins
   * \return The max of val and 0
   */
  IndexType clampLower(IndexType val, const BinSet& AXOM_NOT_USED(binSet)) const
  {
    const IndexType lowerVal = IndexType();
    return val < lowerVal ? lowerVal : val;
  }

private:

  SpatialBoundingBox m_bb;
  LatticeType        m_lattice;
  double             m_expansionFactor;

  ElementSet         m_elementSet;
  BinSet             m_bins[NDIMS];
  BinBitMap          m_binData[NDIMS];

  bool               m_initialized;

};



} // end namespace quest
} // end namespace axom

#endif  // QUEST_IMPLICIT_GRID__HPP_
