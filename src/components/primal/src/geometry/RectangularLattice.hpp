#ifndef PRIMAL_LATTICE__HPP_
#define PRIMAL_LATTICE__HPP_

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
 * \class RectangularLattice
 * \brief A rectangular lattice maps all of space (of dimension NDIMS)
 *        into a rectangular grid of cells identified by integer coordinates.
 *
 * \tparam NDIMS the dimension of the space
 * \tparam SpaceCoordType the type of the spatial coordinates. Default: double
 * \tparam CellCoordType the type of the integral lattice coordinates.
 *         Default: int
 *
 * A rectangular lattice maps every point in space (a SpacePoint) to a unique
 * rectangular cell with integral coordinates (a GridCell). Each lattice cell is
 * identified with the integer grid point at its lowest corner in each dimension.
 * RectangularLattice also maps GridCell coordinates back to spatial coordinates
 * or to the bounding box of the GridCell.
 *
 * A RectangularLattice is defined by an origin (a SpacePoint)
 * and a grid spacing (a SpaceVector).
 */
template < int NDIMS,
           typename SpaceCoordType = double,
           typename CellCoordType = int >
class RectangularLattice
{
public:
  typedef primal::Point< CellCoordType, NDIMS >          GridCell;
  typedef primal::Point< SpaceCoordType, NDIMS >         SpacePoint;
  typedef primal::Vector< SpaceCoordType, NDIMS >        SpaceVector;

  typedef primal::BoundingBox< SpaceCoordType, NDIMS >   SpatialBoundingBox;

  static const SpacePoint s_DefaultOrigin;      /// Default lattice origin
  static const SpaceVector s_DefaultSpacing;    /// Default lattice spacing


  /**
   * \brief Default constructor
   *
   * \note Sets the origin to 0 and spacing to 1 in each dimension
   */
  RectangularLattice()
    : m_origin(SpaceCoordType(0)),
      m_spacing( SpacePoint(SpaceCoordType(1)) ),
      m_invSpacing( SpacePoint(SpaceCoordType(1)) )
  {}

  /**
   * \brief Constructor using a given origin.
   *
   * \param origin The lattice's origin
   * \note Sets the default spacing to 1 in each dimension
   */
  RectangularLattice(const SpacePoint& origin)
    : m_origin(origin),
      m_spacing( SpacePoint(SpaceCoordType(1)) ),
      m_invSpacing( SpacePoint(SpaceCoordType(1)) )
  {}

  /**
   * Constructor using a given origin and spacing along each dimension.
   * \param origin The lattice's origin (Default: s_DefaultOrigin)
   * \param spacing The lattice's spacing (Default: s_DefaultSpacing)
   */
  RectangularLattice(const SpacePoint& origin, const SpaceVector & spacing)
    : m_origin(origin), m_spacing(spacing)
  {
    // Note: m_invSpacing is for efficiency.  It trades divisions for
    // multiplications, and handles dealing with 0-sized spacings
    for (int i=0; i< NDIMS; ++i) {
      m_invSpacing[i] =
       axom::utilities::isNearlyEqual(spacing[i], SpaceCoordType(0))
          ? SpaceCoordType(0)
          : SpaceCoordType(1) / spacing[i];
    }
  }

  /** Accessor for lattice origin   */
  const SpacePoint&  origin() const { return m_origin; }

  /** Accessor for lattice spacing   */
  const SpaceVector& spacing() const { return m_spacing; }

  /** Finds the lattices cell associated with the given space point pt */
  GridCell getGridCell(const SpacePoint& pt) const
  {
    GridCell cell;

    for (int i=0; i< NDIMS; ++i) {
      // Note: Need floor to round down on the negative axes
      cell[i] = static_cast< CellCoordType >(
        std::floor( ( pt[i] - m_origin[i]) * m_invSpacing[i] ) );
    }

    return cell;
  }

  /** Finds the space point associated with the lowest corner of GridCell cell */
  SpacePoint getSpacePoint(const GridCell& cell) const
  {
    SpacePoint pt;
    for (int i=0; i< NDIMS; ++i) {
      pt[i] = m_origin[i] + (m_spacing[i] * cell[i]);
    }

    return pt;
  }

  /**
   * Find the spatial bounding box of a lattice cell
   *
   * \note The bounding box covers all points in space that map to gridCell
   */
  SpatialBoundingBox getCellBounds(const GridCell& cell) const
  {
    return SpatialBoundingBox( getSpacePoint(cell),
                               getSpacePoint(cell.array() +
                                             GridCell(1).array()));
  }

private:
  SpacePoint m_origin;       /// Origin of the grid
  SpaceVector m_spacing;     /// Spacing of the grid
  SpaceVector m_invSpacing;  /// Inverse spacing (for efficiency)
};


/**
 * \brief Helper function to create a Rectangular Lattice from
 * a supplied bounding box and grid resolution
 *
 * \param bbox The bounding box from which to construct a RectangularLattice
 * \param gridRes The resolution of the input grid
 *
 * This is a convenience function to simplify working with bounding boxes.
 * It extracts the lattice spacing from the supplied bounding box
 * and resolution, and sets the lattice origin to the bounding box's
 * minimum corner position.
 */
template < int NDIMS, typename SpaceCoordType, typename CellCoordType >
RectangularLattice< NDIMS, SpaceCoordType, CellCoordType >
rectangular_lattice_from_bounding_box(
  const primal::BoundingBox< SpaceCoordType, NDIMS >& bbox,
  const primal::NumericArray< CellCoordType, NDIMS >& gridRes)
{
  typedef RectangularLattice< NDIMS, SpaceCoordType, CellCoordType > LatticeType;
  typedef typename LatticeType::SpaceVector SpaceVector;

  SpaceVector spacing;

  // Use the resolution and the box range to compute the spacing
  for (int i=0; i< NDIMS; ++i) {
    if (gridRes[i] != SpaceCoordType(0)) {
      spacing[i] = (bbox.getMax()[i] - bbox.getMin()[i]) / gridRes[i];
    }
  }

  return LatticeType(bbox.getMin(), spacing);
}

} // end namespace primal
} // end namespace axom

#endif  // PRIMAL_LATTICE__HPP_
