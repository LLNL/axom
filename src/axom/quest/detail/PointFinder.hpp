/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
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

#ifndef QUEST_POINT_IN_CELL_POINT_FINDER_HPP_
#define QUEST_POINT_IN_CELL_POINT_FINDER_HPP_

#include "axom/quest/geom/ImplicitGrid.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"


namespace axom
{
namespace quest
{

// Predeclare mesh traits class
template <typename mesh_tag>
struct PointInCellTraits;


namespace detail
{

// Predeclare mesh wrapper class
template <typename mesh_tag>
class PointInCellMeshWrapper;



/*!
 * \class PointFinder
 *
 * \brief A class to encapsulate locating points within the cells
 * of a computational mesh.
 *
 * \tparam NDIMS The dimension of the mesh
 * \tparam mesh_tag A tag struct used to identify the mesh
 *
 * \note This class implements part of the functionality of \a PointInCell
 * \note This class assumes the existence of specialized implementations of
 * the following two classes for the provided \a mesh_tag:
 *   \arg axom::quest::PointInCellTraits
 *   \arg axom::quest::detail::PointInCellMeshWrapper
 */
template<int NDIMS, typename mesh_tag>
class PointFinder
{
public:
  typedef quest::ImplicitGrid<NDIMS> GridType;

  typedef typename GridType::SpacePoint SpacePoint;
  typedef typename GridType::SpatialBoundingBox SpatialBoundingBox;

  typedef PointInCellMeshWrapper<mesh_tag> MeshWrapperType;
  typedef typename MeshWrapperType::IndexType IndexType;

public:
  /*!
   * Constructor for PointFinder
   *
   * \param meshWrapper A non-null MeshWrapperType
   * \param res The grid resolution for the spatial acceleration structure
   * \param bboxScaleFactor A number slightly larger than 1 by which to expand
   * cell bounding boxes
   *
   * \sa constructors in PointInCell class for more details about parameters
   */
  PointFinder(const MeshWrapperType* meshWrapper,
              const int* res,
              double bboxScaleFactor)
    : m_meshWrapper(meshWrapper)
  {
    SLIC_ASSERT( m_meshWrapper != nullptr);
    SLIC_ASSERT( bboxScaleFactor >= 1.);

    const int numCells = m_meshWrapper->numElements();

    // setup bounding boxes -- Slightly scaled for robustness

    SpatialBoundingBox meshBBox;
    m_cellBBoxes = std::vector<SpatialBoundingBox>(numCells);
    m_meshWrapper->template computeBoundingBoxes<NDIMS>(bboxScaleFactor,
                                                        m_cellBBoxes, meshBBox);

    // initialize implicit grid, handle case where resolution is a NULL pointer
    if(res != nullptr)
    {
      typedef axom::primal::Point<int, NDIMS> GridResolution;
      GridResolution gridRes(res);
      m_grid.initialize(meshBBox, &gridRes, numCells);
    }
    else
    {
      m_grid.initialize(meshBBox, nullptr, numCells);
    }

    // add mesh elements to grid
    for(int i=0 ; i< numCells ; ++i)
    {
      m_grid.insert( m_cellBBoxes[i], i);
    }
  }

  /*!
   * Query to find the mesh cell containing query point with coordinates \a pos
   *
   * \sa PointInCell::locatePoint() for more details about parameters
   */
  IndexType locatePoint(const double* pos, double* isoparametric) const
  {
    typedef typename GridType::BitsetType BitsetType;

    IndexType containingCell = PointInCellTraits<mesh_tag>::NO_CELL;

    SLIC_ASSERT( pos != nullptr);
    SpacePoint pt(pos);
    SpacePoint isopar;

    // Note: ImplicitGrid::getCandidates() checks the mesh bounding box for us
    BitsetType candidates = m_grid.getCandidates(pt);

    bool foundContainingCell = false;
    for(IndexType cellIdx = candidates.find_first() ;
        !foundContainingCell && cellIdx != BitsetType::npos ;
        cellIdx = candidates.find_next( cellIdx) )
    {
      // First check that pt is in bounding box of element
      if( cellBoundingBox(cellIdx).contains(pt) )
      {
        // if isopar is in the proper range
        if( m_meshWrapper->locatePointInCell(cellIdx, pt.data(),
                                             isopar.data() ) )
        {
          // then we have found the cellID
          foundContainingCell = true;
          containingCell = cellIdx;
        }
      }
    }

    // Copy data back to input parameter isoparametric, if necessary
    if(isoparametric != nullptr)
    {
      isopar.array().to_array(isoparametric);
    }

    return containingCell;
  }

  /*! Returns a const reference to the given cells's bounding box */
  const SpatialBoundingBox& cellBoundingBox(IndexType cellIdx) const
  {
    return m_cellBBoxes[cellIdx];
  }


private:
  GridType m_grid;
  const MeshWrapperType* m_meshWrapper;
  std::vector<SpatialBoundingBox> m_cellBBoxes;
};


} // end namespace detail
} // end namespace quest
} // end namespace axom


#endif // QUEST_POINT_IN_CELL_POINT_FINDER_HPP_
