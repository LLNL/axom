#ifndef PRIMAL_GRID__HPP_
#define PRIMAL_GRID__HPP_

#include "axom/config.hpp"

#include "fmt/fmt.hpp"
#include "slic/slic.hpp"

#include "primal/BoundingBox.hpp"
#include "primal/Point.hpp"
#include "primal/Vector.hpp"

#include <cmath>  // for std::floor

namespace axom {
namespace primal {

/**
 * \class Grid
 * \brief A grid quantizes points in space to integral grid cells based on an origin
 *        and grid spacing.
 *
 * Grids are unbounded quotient spaces that map every point in space (a SpacePoint)
 * to a single GridCell with integral coordinates. Grids are aligned with the
 * Euclidean axes and are defined by an origin (SpacePoint) and a grid spacing
 * (SpaceVector). Grid cells are identified with the grid point at the lowest corner
 * in each dimension.
 */
template < int NDIMS,
           typename SpaceCoordType = double,
           typename GridCoordType = int >
class Grid
{
public:
  typedef primal::Point< GridCoordType, NDIMS >          GridCell;
  typedef primal::Point< SpaceCoordType, NDIMS >         SpacePoint;
  typedef primal::Vector< SpaceCoordType, NDIMS >        SpaceVector;

  typedef primal::BoundingBox< SpaceCoordType, NDIMS >   SpatialBoundingBox;

  static const SpacePoint s_DefaultOrigin;      /// Default origin for a grid
  static const SpaceVector s_DefaultSpacing;    /// Default spacing for a grid

  /**
   * Constructor for a Grid with a given origin and spacing along each dimension.
   */
  Grid(const SpacePoint& origin = s_DefaultOrigin,
       const SpaceVector & spacing = s_DefaultSpacing ):
    m_origin(origin), m_spacing(spacing)
  {
    for (int i=0; i< NDIMS; ++i) {
      SLIC_ASSERT( spacing[i] != SpaceCoordType(0) );
      m_invSpacing[i] = SpaceCoordType(1) / spacing[i];
    }
  }

  /** Accessor for grid origin   */
  const SpacePoint&  gridOrigin() const { return m_origin; }

  /** Accessor for grid spacing   */
  const SpaceVector& gridSpacing() const { return m_spacing; }

  /** Finds the grid cell associated with the given space point pt */
  GridCell getGridCell(const SpacePoint& pt) const
  {
    GridCell gridCell;

    for (int i=0; i< NDIMS; ++i) {
      // Note: Need floor to round down on the negative axes
      gridCell[i] = static_cast< GridCoordType >(
        std::floor( ( pt[i] - m_origin[i]) * m_invSpacing[i] ) );
    }

    return gridCell;
  }

  /** Finds the space point associated with the lowest corner of GridCell cell */
  SpacePoint getSpacePoint(const GridCell& cell) const
  {
    SpacePoint pt;
    for (int i=0; i< NDIMS; ++i) {
      pt[i] = m_origin[i] + m_spacing[i] * cell[i];
    }

    return pt;
  }

  /**
   * Find the spatial bounding box of a grid cellcell
   *
   * \note The bounding box covers all points in space that map to gridCell
   */
  SpatialBoundingBox getCellBounds(const GridCell& gridCell) const
  {
    return SpatialBoundingBox( getSpacePoint(gridCell),
                               getSpacePoint(gridCell.array() +
                                             GridCell(1).array()));
  }

private:
  SpacePoint m_origin;       /// Origin of the grid
  SpaceVector m_spacing;     /// Spacing of the grid
  SpaceVector m_invSpacing;  /// Inverse spacing for efficiency
};

// Initialize the default origin to zeros in each dimension
template < int NDIMS, typename SpaceCoordType, typename GridCoordType >
const typename Grid< NDIMS, SpaceCoordType, GridCoordType >::SpacePoint
Grid< NDIMS, SpaceCoordType,GridCoordType >::s_DefaultOrigin
  = SpacePoint::zero();

// Initialize the default spacing to ones in each dimension
template < int NDIMS, typename SpaceCoordType, typename GridCoordType >
const typename Grid< NDIMS, SpaceCoordType, GridCoordType >::SpaceVector
Grid< NDIMS, SpaceCoordType,GridCoordType >::s_DefaultSpacing
  = SpaceVector( SpacePoint(SpaceCoordType(1)));

/**
 * Helper function to screate a Grid from a bounding box and a grid resolution
 * \note This is a convenience function to simplify working with bounding boxes
 */
template < int NDIMS, typename SpaceCoordType, typename GridCoordType >
Grid< NDIMS, SpaceCoordType, GridCoordType > make_grid_from_bounding_box(
  const primal::BoundingBox< SpaceCoordType, NDIMS >& bbox,
  const primal::NumericArray< GridCoordType, NDIMS >& gridRes)
{
  typedef Grid< NDIMS, SpaceCoordType, GridCoordType > GridType;
  typedef typename GridType::SpaceVector SpaceVector;

  SpaceVector spacing;
  for (int i=0; i< NDIMS; ++i) {
    if (gridRes[i] != SpaceCoordType(0)) {
      spacing[i] = (bbox.getMax()[i] - bbox.getMin()[i]) / gridRes[i];
    }
  }

  return GridType(bbox.getMin(), spacing);
}

} // end namespace primal
} // end namespace axom

#endif  // PRIMAL_GRID__HPP_
