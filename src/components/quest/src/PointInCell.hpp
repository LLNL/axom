/*
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


#ifndef QUEST_POINT_IN_CELL_HPP_
#define QUEST_POINT_IN_CELL_HPP_


#include "axom/config.hpp"
#include "axom/Macros.hpp"

#include "primal/Point.hpp"
#include "primal/BoundingBox.hpp"

#include "quest/ImplicitGrid.hpp"

#ifdef AXOM_USE_MFEM
# include "quest/PointInCell_impl_mfem.hpp"
#endif

namespace axom {
namespace quest {

namespace detail {

/*!
 * \class PointInCellTraits
 * \brief A traits class for the mesh associated with a PointInCell query.
 *
 * \tparam mesh_tag A tag struct used to identify the mesh
 *
 * An implementation of \a PointInCellTraits on \a mesh_tag must define:
 * \arg A \a MeshType (e.g. mfem::Mesh or axom::mint::Mesh)
 * \arg An \a IndexType (e.g. a signed or unsigned integral type)
 * \arg An \a IndexType variable named \a NO_CELL (e.g. with value -1 for signed ints),
 * indicating an invalid cell index in the mesh
 */
template<typename mesh_tag>
struct PointInCellTraits;

/*!
 * \class PointInCellMeshWrapper
 * \brief A wrapper class for accessing the mesh instance
 * associated with a PointInCell query
 *
 * \tparam mesh_tag A tag struct used to identify the mesh
 *
 * An implementation of a \a PointInCellMeshWrapper on \a mesh_tag
 * must expose the following API:
 * \arg template<DIM>
 *      void computeBoundingBox(
 *          double, const std::vector< BoundingBox<double, DIM>&,
 *          const BoundingBox<double, DIM>  &) const;
 * \arg void reconstructPoint(IndexType, const double*, double*) const;
 * \arg bool locatePointInCell(IndexType, const double*, double*) const;
 * \arg int numElements() const;
 * \arg int meshDimension() const;
 */
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
 *   \arg axom::quest::detail::PointInCellTraits
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
   *
   * \sa constructors in PointInCell class for more details about parameters
   */
  PointFinder(MeshWrapperType* meshWrapper, int* res)
    : m_meshWrapper(meshWrapper)
  {
    SLIC_ASSERT( m_meshWrapper != AXOM_NULLPTR);

    const int numCells = m_meshWrapper->numElements();

    // setup bounding boxes -- Slightly scaled for robustness

    SpatialBoundingBox meshBBox;
    m_cellBBoxes = std::vector<SpatialBoundingBox>(numCells);
    double EPS = 1e-8;
    double bboxScaleFactor = 1 + EPS;
    m_meshWrapper->template computeBoundingBoxes<NDIMS>(bboxScaleFactor, m_cellBBoxes, meshBBox);

    // initialize implicit grid
    typedef axom::primal::Point<int, NDIMS> GridResolution;
    GridResolution gridRes;
    if(res != AXOM_NULLPTR)
    {
      gridRes = GridResolution(res);
    }
    m_grid.initialize(meshBBox, gridRes, numCells);

    // add mesh elements to grid
    for(int i=0; i< numCells; ++i)
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

    SLIC_ASSERT( pos != AXOM_NULLPTR);
    SpacePoint pt(pos);
    SpacePoint isopar;

    // Note: ImplicitGrid::getCandidates() checks the mesh bounding box for us
    BitsetType candidates = m_grid.getCandidates(pt);

    bool foundContainingCell = false;
    for(std::size_t cellIdx = candidates.find_first();
        !foundContainingCell && cellIdx != BitsetType::npos;
        cellIdx = candidates.find_next( cellIdx) )
    {
      // First check that pt is in bounding box of element
      if( cellBoundingBox(cellIdx).contains(pt) )
      {
        // if isopar is in the proper range
        if( m_meshWrapper->locatePointInCell(cellIdx, pt.data(), isopar.data() ) )
        {
          // then we have found the cellID
          foundContainingCell = true;
          containingCell = cellIdx;
        }
      }
    }

    // Copy data back to input parameter isoparametric, if necessary
    if(isoparametric != AXOM_NULLPTR)
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
  GridType    m_grid;
  MeshWrapperType* m_meshWrapper;
  std::vector<SpatialBoundingBox> m_cellBBoxes;
};

} // end namespace detail


/*!
 * \class PointInCell
 *
 * \brief A class to accelerate Point-In-Cell queries on a computational mesh.
 *
 * A point in cell query over a computational mesh determines if a given point is
 * contained in one of its cells.  If so, it returns the index of the cell
 * containing the point as well as the isoparametric coordinates of the point
 * within the cell.
 *
 * \tparam mesh_tag A tag type for the the underlying computational mesh
 * There must also be corresponding template specializations of
 * \a PointInCellMeshWrapper<mesh_tag> and \a PointInCellTraits<mesh_tag> for
 * provided mesh_tag \a axom::quest::detail namespace
 *
 * \note This class was designed to support point in cell queries against
 * 2D or 3D computational meshes. The queries to the mesh are wrapped in a
 * PointInCellMeshWrapper class templated on a \a mesh_tag type.
 * To extend this design, one must create a new mesh tag (e.g.. an empty struct
 * \a custom_mesh_tag) and add custom template specializations of
 * PointInCellMeshWrapper and PointInCellTraits for this tag in the
 * axom::quest::detail namespace.
 *
 * \sa PointInCell_impl_mfem.hpp for a specialized implementation
 * for [mfem](http://mfem.org) meshes of arbitrary order.
 */
template<typename mesh_tag>
class PointInCell
{
public:
  typedef detail::PointInCellTraits<mesh_tag> MeshTraits;
  typedef typename MeshTraits::MeshType MeshType;
  typedef typename MeshTraits::IndexType IndexType;

  typedef detail::PointInCellMeshWrapper<mesh_tag> MeshWrapperType;

  typedef detail::PointFinder<2, mesh_tag> PointFinder2D;
  typedef detail::PointFinder<3, mesh_tag> PointFinder3D;


  /*!
   * Construct a point in cell query structure over a computational mesh
   *
   * \param mesh A pointer to the computational mesh
   * \param resolution Grid resolution for the spatial index.
   *
   * \note If the resolution is not provided, a heuristic based on the number
   * of cells in the mesh is used to set the resolution.
   * \sa ImplicitGrid
   *
   * \pre \a mesh must not be a NULL pointer
   * \pre If resolution is not NULL, it must have space for at least
   * meshDimension() entries.
   */
  PointInCell(MeshType* mesh, int* resolution = AXOM_NULLPTR)
    : m_meshWrapper(mesh), m_pointFinder2D(AXOM_NULLPTR), m_pointFinder3D(AXOM_NULLPTR)
  {
    SLIC_ASSERT(mesh != AXOM_NULLPTR);

    // Allocate a 2D or 3D PointFinder instance, depending on mesh dimension.
    // Note: Only one of these will be allocated in a PointInCell instance
    switch(m_meshWrapper.meshDimension())
    {
    case 2:
      m_pointFinder2D = new PointFinder2D(&m_meshWrapper, resolution);
      break;
    case 3:
      m_pointFinder3D = new PointFinder3D(&m_meshWrapper, resolution);
      break;
    default:
      SLIC_ERROR("Point in Cell query only defined for 2D or 3D meshes.");
      break;
    }
  }

  ~PointInCell()
  {
    if( m_pointFinder2D != AXOM_NULLPTR)
    {
      delete m_pointFinder2D;
      m_pointFinder2D = AXOM_NULLPTR;
    }

    if( m_pointFinder3D != AXOM_NULLPTR)
    {
      delete m_pointFinder3D;
      m_pointFinder3D = AXOM_NULLPTR;
    }
  }

  /*!
   * Attempt to find the index of the mesh cell containing the given point.
   *
   * If found, and \a isopar is not NULL, \a isopar contains the isoparametric
   * coordinates the point within this cell.
   *
   * \param[in] pos The coordinates of the query point in space
   * \param[out] isopar The isoparametric coordinates of the query pt.
   * Only valid when a cell is found.
   *
   * \return The index of the mesh cell containing the query point.
   * If no cell is found, returns the special value \a MeshTraits::NO_CELL
   *
   * \pre \a pos is a non-null array with at least \a meshDimension() coordinates
   * \pre When not NULL, \a isopar has space for at least \a meshDimension() coordinates
   */
  IndexType locatePoint(const double* pos, double* isopar = AXOM_NULLPTR) const
  {
    SLIC_ASSERT(pos != AXOM_NULLPTR);

    IndexType cellIndex = MeshTraits::NO_CELL;

    switch(m_meshWrapper.meshDimension())
    {
    case 2:
      cellIndex = m_pointFinder2D->locatePoint(pos, isopar);
      break;
    case 3:
      cellIndex = m_pointFinder3D->locatePoint(pos, isopar);
      break;
    default:
      SLIC_ERROR("Point in Cell query only defined for 2D or 3D meshes.");
      break;
    }

    return cellIndex;
  }

  /*!
   *  Determine if a query point is located within a specified mesh cell
   *
   *  \param [in] cellIdx The index of the cell within the mesh
   *  \param [in] pos The coordinates of the query point in space
   *  \param [out] isopar The isoparametric coordinates of the point
   *  within cell \a cellIdx.  Only valid when the return value is true
   *
   *  \return True, if the query point is located within the specified cell
   *  \pre \a pos is not NULL and has \a meshDimension() entries
   *  \pre \a isopar is not NULL and has space for \a meshDimension() entries
   */
  bool locatePointInCell(IndexType cellIdx, const double* pos, double* isopar) const
  {
    // Early return if point is not within cell's bounding box
    if(! withinBoundingBox(cellIdx, pos) )
    {
      return false;
    }

    return m_meshWrapper.locatePointInCell(cellIdx, pos, isopar);
  }


  /*!
   * Evaluate the position of a point within a mesh cell at the given
   * isoparametric coordinates.
   *
   * \param [in] cellIdx The index of the cell within the mesh
   * \param [in[ isopar The isoparametric coordinates at which to evaluate
   * \param [out] pos The computed coordinates of the evaluated point
   */
  void reconstructPoint(IndexType cellIdx, const double* isopar, double* pos) const
  {
    m_meshWrapper.reconstructPoint(cellIdx, isopar, pos);
  }

  /*! Returns the dimension of the mesh */
  int meshDimension() const { return m_meshWrapper.meshDimension(); }

private:

  /*!
   * Utility function to check the given point against an element's bounding box
   * \param [in] cellIdx Index of the cell within the mesh
   * \param [in] pos Position of the point in space
   *
   * \return True if the point is contained in the cell's bounding box
   */
  bool withinBoundingBox(IndexType cellIdx, const double* pos) const
  {
    typedef axom::primal::Point<double, 2> Point2D;
    typedef axom::primal::Point<double, 3> Point3D;

    switch(meshDimension())
    {
    case 2:
      return m_pointFinder2D->elementBoundingBox(cellIdx).contains(Point2D(pos));
    case 3:
      return m_pointFinder3D->elementBoundingBox(cellIdx).contains(Point3D(pos));
    }
    return false;
  }


private:
  MeshWrapperType m_meshWrapper;

  PointFinder2D* m_pointFinder2D;
  PointFinder3D* m_pointFinder3D;
  
  // double m_eps; // for scaling the bounding boxes? 
};




} // end namespace quest
} // end namespace axom

#endif // QUEST_POINT_IN_CELL_HPP_
