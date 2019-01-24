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

#ifndef PRIMAL_RECTANGULAR_LATTICE_HPP_
#define PRIMAL_RECTANGULAR_LATTICE_HPP_

#include "axom/config.hpp"
#include "axom/core/Types.hpp"

#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"

#include <cmath>     // for std::floor
#include <iostream>  // for ostream

namespace axom
{
namespace primal
{

/*!
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
 * identified with the integer grid point at its lowest corner in each
 * dimension.
 *
 * RectangularLattice also maps GridCell coordinates back to spatial coordinates
 * or to the bounding box of the GridCell.
 *
 * GridCells follow a half-open boundary convention. Points on their lower
 * boundaries get mapped to the cell, while points on their upper boundaries
 * get mapped to neighboring cells.
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

  /*!
   * \brief Default constructor
   *
   * \note Sets the origin to 0 and spacing to 1 in each dimension
   */
  RectangularLattice()
    : m_origin(SpaceCoordType(0)),
    m_spacing( SpacePoint(SpaceCoordType(1)) ),
    m_invSpacing( SpacePoint(SpaceCoordType(1)) )
  {}

  /*!
   * \brief Constructor from a given origin.
   *
   * \param origin The lattice's origin
   * \note Spacing will be set to default (vector of ones)
   */
  RectangularLattice(const SpacePoint& origin)
    : m_origin(origin),
    m_spacing( SpacePoint(SpaceCoordType(1)) ),
    m_invSpacing( SpacePoint(SpaceCoordType(1)) )
  {}

  /*!
   * \brief Constructor from a given origin and spacing.
   *
   * \param origin The lattice's origin
   * \param spacing The lattice's spacing
   */
  RectangularLattice(const SpacePoint& origin, const SpaceVector & spacing)
    : m_origin(origin), m_spacing(spacing)
  {
    // Note: m_invSpacing is for efficiency.  It trades divisions for
    // multiplications, and handles dealing with 0-sized spacings
    for (int i=0 ; i< NDIMS ; ++i)
    {
      m_invSpacing[i] =
        axom::utilities::isNearlyEqual(spacing[i], SpaceCoordType(0))
        ? SpaceCoordType(0)
        : SpaceCoordType(1) / spacing[i];
    }
  }

  /*!
   * \brief Constructor from SpaceCoordType arrays
   *
   * \param origin_date An array containing the origin's coordinates
   * \param spacing_data An array containing the spacing coordinates
   *
   * \pre When origin_data is not NULL, it must have at least NDIMS entries
   * \note Origin will be set to the zero vector if origin_data pointer is NULL
   *
   * \pre When spacing_data is not NULL, it must have at least NDIMS entries
   * \note Spacing will be set to vector or ones if pointer is NULL
   */
  RectangularLattice(SpaceCoordType* origin_data,
                     SpaceCoordType* spacing_data)
  {
    m_origin = (origin_data != nullptr)
               ? SpacePoint(origin_data)
               : SpacePoint::zero();

    m_spacing = (spacing_data != nullptr)
                ? SpaceVector(spacing_data)
                : SpaceVector(SpaceCoordType(1));

    // Note: m_invSpacing is for efficiency.  It trades divisions for
    // multiplications, and handles dealing with 0-sized spacings
    for (int i=0 ; i< NDIMS ; ++i)
    {
      m_invSpacing[i] =
        axom::utilities::isNearlyEqual(m_spacing[i], SpaceCoordType(0))
        ? SpaceCoordType(0)
        : SpaceCoordType(1) / m_spacing[i];
    }
  }

  /*! Accessor for lattice origin   */
  const SpacePoint&  origin() const { return m_origin; }

  /*! Accessor for lattice spacing   */
  const SpaceVector& spacing() const { return m_spacing; }

  /*! Returns the lattice cell associated with the given space point pt */
  GridCell gridCell(const SpacePoint& pt) const
  {
    GridCell cell;

    for (int i=0 ; i< NDIMS ; ++i)
    {
      // Note: Always round down to negative infinity
      cell[i] = static_cast< CellCoordType >(
        std::floor( (pt[i] - m_origin[i]) * m_invSpacing[i] ) );
    }

    return cell;
  }

  /*! Returns the space point associated with lowest corner of GridCell cell */
  SpacePoint spacePoint(const GridCell& cell) const
  {
    SpacePoint pt;
    for (int i=0 ; i< NDIMS ; ++i)
    {
      pt[i] = m_origin[i] + (m_spacing[i] * cell[i]);
    }

    return pt;
  }

  /*!
   * \brief Find the spatial bounding box of a lattice cell
   *
   * \note The bounding box covers all points in space that map to gridCell
   */
  SpatialBoundingBox cellBounds(const GridCell& cell) const
  {
    return SpatialBoundingBox( spacePoint(cell),
                               spacePoint(cell.array() +
                                          GridCell(1).array()));
  }

  /*! Simple formatted print of a rectangular lattice */
  std::ostream& print(std::ostream& os) const
  {
    os <<"{ origin:"<<m_origin
       <<"; spacing:"<< m_spacing <<" }";
    return os;
  }

private:
  SpacePoint m_origin;       /// Origin of the grid
  SpaceVector m_spacing;     /// Spacing of the grid
  SpaceVector m_invSpacing;  /// Inverse spacing (for efficiency)
};

/*!
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
  typedef RectangularLattice< NDIMS, SpaceCoordType,
                              CellCoordType > LatticeType;
  typedef typename LatticeType::SpaceVector SpaceVector;

  SpaceVector spacing;

  // Use the resolution and the box range to compute the spacing
  for (int i=0 ; i< NDIMS ; ++i)
  {
    if (gridRes[i] != SpaceCoordType(0))
    {
      spacing[i] = (bbox.getMax()[i] - bbox.getMin()[i]) / gridRes[i];
    }
  }

  return LatticeType(bbox.getMin(), spacing);
}

/// --------------------------------------------------------------------------------
/// Free functions implementing comparison and arithmetic operators
/// --------------------------------------------------------------------------------

/*! Equality operator on two RectangularLattices */
template < int NDIMS, typename SpaceCoordType, typename CellCoordType >
bool operator==(
  const RectangularLattice< NDIMS,SpaceCoordType,CellCoordType >& lhs,
  const RectangularLattice< NDIMS,SpaceCoordType,CellCoordType >& rhs )
{
  return lhs.origin() == rhs.origin() && lhs.spacing() == rhs.spacing();
}

/*! Inequality operator on two RectangularLattices */
template < int NDIMS, typename SpaceCoordType, typename CellCoordType >
bool operator!=(
  const RectangularLattice< NDIMS,SpaceCoordType,CellCoordType >& lhs,
  const RectangularLattice< NDIMS,SpaceCoordType,CellCoordType >& rhs )
{
  return !(lhs == rhs);
}

/*! Stream output operator on a RectangularLattice */
template < int NDIMS, typename SpaceCoordType, typename CellCoordType >
std::ostream& operator<<(
  std::ostream & os,
  const RectangularLattice< NDIMS, SpaceCoordType, CellCoordType > & lattice)
{
  return lattice.print(os);
}

} // end namespace primal
} // end namespace axom

#endif  // PRIMAL_RECTANGULAR_LATTICE_HPP_
