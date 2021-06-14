// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_UNIFORMGRID_HPP_
#define AXOM_SPIN_UNIFORMGRID_HPP_

#include "axom/core/utilities/Utilities.hpp"

#include "axom/slic/interface/slic.hpp"

#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/NumericArray.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/spin/RectangularLattice.hpp"

// C/C++ includes
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>
#include <stdio.h>

namespace axom
{
namespace spin
{
/*!
 * \class UniformGrid
 *
 * \brief A spatial index defined by origin, spacing, and resolution.
 *
 * \tparam T The type of object that this UniformGrid will hold
 * \tparam NDIMS The UniformGrid dimensionality: supported values are 2 or 3
 *
 * The UniformGrid class is parameterized on the type of object it will contain
 * and the dimensionality of the grid.  Currently the UniformGrid will index
 * 2D and 3D space.
 *
 * The UniformGrid divides space into a number of bins, extending in a block
 * from an origin point and arranged in row-major order.  Each bin extends over
 * an axis-aligned rectangle or parallelepiped.  The size of a bin is specified
 * by the spacing constructor argument, and the number of bins is specified by
 * the res constructor argument.
 *
 * To use this class, a user will instantiate a UniformGrid and insert a number
 * of objects (each with its own bounding box).  The insert operation puts the
 * object into some of the bins.  The user may retrieve the bin index of any
 * point, then retrieve any objects associated with the bin at that index.
 */
template <typename T, int NDIMS>
class UniformGrid
{
public:
  /*! \brief The type used for specifying spatial extent of the contents */
  using BoxType = primal::BoundingBox<double, NDIMS>;

  /*! \brief The type used to query the index */
  using PointType = primal::Point<double, NDIMS>;

private:
  /*! \brief The type used for mapping points in space to grid cells */
  using LatticeType = RectangularLattice<NDIMS, double, int>;

  /*! \brief The type used to represent integer-valued grid cells */
  using GridCell = typename LatticeType::GridCell;

public:
  /*!
   * \brief Constructor specifying bounding box (min and max points) and
   *        number of bins.
   *
   * The lower and upper bounds are specified as double arrays.  Each argument
   * points to an array of at least length NDIMS.  The nth element in each array
   * specifies, for the nth dimension, the lower bound (inclusive), the upper
   * bound (exclusive), and the number of bins.  Thus,
   *     UniformGrid({0., 0.}, {2., 4.}, {2, 4});
   * defines an index extending from 0 up to 2 divided into two bins in the
   * x-dimension and from 0 up to 4 divided into four bins in the y dimension.
   */
  UniformGrid(const double* lower_bound, const double* upper_bound, const int* res);

  /*!
   * \brief Constructor specifying bounding box and number of bins.
   */
  UniformGrid(const BoxType& bbox, const int* res);

  /*! \brief Destructor: present for symmetry with constructor */
  ~UniformGrid();

  /*!
   * \brief Returns the index of the bin containing the specified point.
   *
   * If the specified point is outside the grid, the function returns a
   * special value (INVALID_BIN_INDEX).  Note that the upper boundaries of a
   * UniformGrid are not in a bin.  Thus, if a UniformGrid divides the region
   * from (0., 0.) to (10., 10.) into 10 x 10 bins, the code
   *     getBinIndex(PointType::makePoint(10., 10.))
   * will return INVALID_BIN_INDEX,
   * \param [in] pt The point to query.
   */
  int getBinIndex(const PointType& pt) const;

  /*! \brief Returns the number of bins in this UniformGrid. */
  int getNumBins() const;

  /*!
   * \brief Returns whether the bin specified by index is empty.
   *
   * Returns true if index is invalid.
   * \param [in] index The index of the bin to test.
   */
  bool isBinEmpty(int index) const;

  /*!
   * \brief Returns the contents of the bin indicated by index.
   *
   * It is an error if index is invalid.
   * \param [in] index The index of the bin to retrieve.
   */
  std::vector<T>& getBinContents(int index);

  /*!
   * \brief Returns the contents of the bin indicated by index.
   *
   * It is an error if index is invalid.  This is the const version.
   * \param [in] index The index of the bin to retrieve.
   */
  const std::vector<T>& getBinContents(int index) const;

  /*!
   * \brief Returns true if index is valid; that is, refers to a valid bin.
   *
   * The region of space indexed by this UniformGrid is divided into bins.
   * Each bin has a unique integer index.  If an integer is greater than or
   * equal to zero and less than or equal to the greatest bin index, this
   * method will return true.  Otherwise, the integer does not refer to any
   * bin and this method will return false.
   */
  bool isValidIndex(int index) const;

  /*!
   * \brief Returns the indices of the bins for a given bounding box.
   *
   * This method returns the indices of the bins that intersect BB.  Any part
   * of BB that falls outside this UniformGrid is disregarded.  If BB falls
   * completely outside this UniformGrid, this method returns an empty list.
   */
  const std::vector<int> getBinsForBbox(const BoxType& BB) const;

  /*!
   * \brief Clears the bin indicated by index.
   *
   * No-op if index is invalid.
   * \param [in] index The index of the bin to clear.
   */
  void clear(int index);

  /*!
   * \brief Inserts obj into each bin overlapped by BB.
   *
   * No error is signalled if BB falls partly or wholly outside the UniformGrid.
   *
   * \param [in] BB The region in which to record obj
   * \param [in] obj The object to insert into any bins overlapped by BB
   */
  void insert(const BoxType& BB, const T& obj);

  /*!
   * \brief A special value indicating any location not in the UniformGrid.
   *
   * Returned by getBinIndex() if that function is passed an exterior point.
   * An error is thrown if getBinContents() is passed INVALID_BIN_INDEX.
   * clear(INVALID_BIN_INDEX) is a no-op; isBinEmpty(INVALID_BIN_INDEX) returns
   * true.
   */
  enum
  {
    INVALID_BIN_INDEX = -1
  };

private:
  /*!
   *****************************************************************************
   * \brief Returns the closest UniformGrid cell to the query point pt.
   *
   * When pt lies within an integer valued grid cell, that cell is returned.
   * When pt is on the boundary, or outside the bounding box of the
   * UniformGrid, the returned grid cell is clamped to its closest grid cell.
   * I.e. each coordinate i of the returned GridCell is between
   * zero and (m_resolution[i] -1).
   *
   * \param [in] pt The query point
   * \return The integer valued grid cell closest to pt
   *****************************************************************************
   */
  GridCell getClampedGridCell(const PointType& pt) const;

  /*! \brief Adds an object obj to the bin at index index */
  void addObj(const T& obj, int index);

  /*!
   *****************************************************************************
   * \brief Common constructor code (until we can use delegating ctors)
   *
   * res specifies the resolution in each dimension.  If nullptr,
   * default_res is used in each dimension.
   *****************************************************************************
   */
  void initialize(const int* res, int default_res = 1);

protected:
  /*! \brief The default constructor should not be used, so it is protected. */
  UniformGrid();

private:
  BoxType m_boundingBox;
  LatticeType m_lattice;

  int m_resolution[NDIMS];
  int m_strides[NDIMS];

  struct Bin
  {
    std::vector<T> ObjectArray;
    //other stuff
  };
  std::vector<Bin> m_bins;

  DISABLE_COPY_AND_ASSIGNMENT(UniformGrid);
  DISABLE_MOVE_AND_ASSIGNMENT(UniformGrid);

};  //end class

}  //end namespace spin
}  //end namespace axom

//------------------------------------------------------------------------------
// UniformGrid implementation
//------------------------------------------------------------------------------
namespace axom
{
namespace spin
{
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
UniformGrid<T, NDIMS>::UniformGrid()
  : m_boundingBox(PointType::zero(), PointType(100))
  , m_lattice(PointType::zero(), PointType(1))
{
  SLIC_ASSERT((NDIMS == 3) || (NDIMS == 2));

  const int DEFAULT_RES = 100;

  initialize(nullptr, DEFAULT_RES);
}

template <typename T, int NDIMS>
UniformGrid<T, NDIMS>::UniformGrid(const double* lower_bound,
                                   const double* upper_bound,
                                   const int* res)
{
  SLIC_ASSERT(lower_bound != nullptr);
  SLIC_ASSERT(upper_bound != nullptr);
  SLIC_ASSERT(res != nullptr);
  SLIC_ASSERT((NDIMS == 3) || (NDIMS == 2));

  // set up the bounding box for point conversions
  m_boundingBox = BoxType(PointType(lower_bound), PointType(upper_bound));

  initialize(res);

  // set up the lattice for point conversions
  m_lattice = rectangular_lattice_from_bounding_box(
    m_boundingBox,
    primal::NumericArray<int, NDIMS>(m_resolution));
}

template <typename T, int NDIMS>
UniformGrid<T, NDIMS>::UniformGrid(const BoxType& bbox, const int* res)
  : m_boundingBox(bbox)
{
  SLIC_ASSERT(res != nullptr);
  SLIC_ASSERT((NDIMS == 3) || (NDIMS == 2));

  initialize(res);

  // set up the bounding box and lattice for point conversions
  m_lattice = rectangular_lattice_from_bounding_box(
    m_boundingBox,
    primal::NumericArray<int, NDIMS>(m_resolution));
}

template <typename T, int NDIMS>
void UniformGrid<T, NDIMS>::initialize(const int* res, int default_res)
{
  // set up the grid resolution and (row-major) strides
  m_resolution[0] = (res ? res[0] : default_res);
  m_strides[0] = 1;
  for(int i = 1; i < NDIMS; ++i)
  {
    m_resolution[i] = (res ? res[i] : default_res);
    m_strides[i] = m_strides[i - 1] * m_resolution[i - 1];
  }

  // initialize space for the bins
  const int numBins = m_strides[NDIMS - 1] * m_resolution[NDIMS - 1];
  m_bins.resize(numBins);

  // scale the bounding box by a little to account for boundaries
  const double EPS = 1e-12;
  m_boundingBox.scale(1. + EPS);
}

template <typename T, int NDIMS>
UniformGrid<T, NDIMS>::~UniformGrid()
{ }

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
int UniformGrid<T, NDIMS>::getBinIndex(const PointType& pt) const
{
  // Index is only valid when the point is within the bounding box
  if(!m_boundingBox.contains(pt))
  {
    return INVALID_BIN_INDEX;
  }

  // Find pt's integer grid cell within the uniform grid
  GridCell cell = getClampedGridCell(pt);

  // compute the linear index (row-major) of the cell
  int res = cell[0];
  for(int i = 1; i < NDIMS; ++i)
  {
    res += m_strides[i] * cell[i];
  }

  return res;
}

template <typename T, int NDIMS>
bool UniformGrid<T, NDIMS>::isValidIndex(int index) const
{
  return index >= 0 && index < static_cast<int>(m_bins.size());
}

template <typename T, int NDIMS>
int UniformGrid<T, NDIMS>::getNumBins() const
{
  return static_cast<int>(m_bins.size());
}

template <typename T, int NDIMS>
bool UniformGrid<T, NDIMS>::isBinEmpty(int index) const
{
  if(!isValidIndex(index))
  {
    return true;
  }

  return m_bins[index].ObjectArray.empty();
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::vector<T>& UniformGrid<T, NDIMS>::getBinContents(int index)
{
  SLIC_ASSERT(isValidIndex(index));

  return m_bins[index].ObjectArray;
}

template <typename T, int NDIMS>
const std::vector<T>& UniformGrid<T, NDIMS>::getBinContents(int index) const
{
  SLIC_ASSERT(isValidIndex(index));

  return m_bins[index].ObjectArray;
}

template <typename T, int NDIMS>
const std::vector<int> UniformGrid<T, NDIMS>::getBinsForBbox(const BoxType& BB) const
{
  std::vector<int> retval;

  if(!m_boundingBox.intersectsWith(BB))
  {
    return retval;
  }

  const GridCell lowerCell = getClampedGridCell(BB.getMin());
  const GridCell upperCell = getClampedGridCell(BB.getMax());

  // Recall that NDIMS is 2 or 3
  const int kLower = (NDIMS == 2) ? 0 : lowerCell[2];
  const int kUpper = (NDIMS == 2) ? 0 : upperCell[2];
  const int kStride = (NDIMS == 2) ? 1 : m_strides[2];

  for(int k = kLower; k <= kUpper; ++k)
  {
    const int kOffset = k * kStride;
    for(int j = lowerCell[1]; j <= upperCell[1]; ++j)
    {
      const int jOffset = j * m_strides[1] + kOffset;
      for(int i = lowerCell[0]; i <= upperCell[0]; ++i)
      {
        retval.push_back(i + jOffset);
      }
    }
  }

  return retval;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
void UniformGrid<T, NDIMS>::clear(int index)
{
  if(isValidIndex(index))
  {
    m_bins[index].ObjectArray.clear();
  }
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
void UniformGrid<T, NDIMS>::addObj(const T& obj, int index)
{
  SLIC_CHECK(isValidIndex(index));

  if(isValidIndex(index))
  {
    m_bins[index].ObjectArray.push_back(obj);
  }
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
void UniformGrid<T, NDIMS>::insert(const BoxType& BB, const T& obj)
{
  SLIC_ASSERT((NDIMS == 3) || (NDIMS == 2));

  const std::vector<int> bidxs = getBinsForBbox(BB);

  const int numBins = static_cast<int>(bidxs.size());
  for(int i = 0; i < numBins; ++i)
  {
    addObj(obj, bidxs[i]);
  }
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
typename UniformGrid<T, NDIMS>::GridCell UniformGrid<T, NDIMS>::getClampedGridCell(
  const PointType& pt) const
{
  GridCell cell = m_lattice.gridCell(pt);

  // clamp the grid coordinates to lie in a cell of the grid
  for(int i = 0; i < NDIMS; ++i)
  {
    cell[i] = axom::utilities::clampVal(cell[i], 0, m_resolution[i] - 1);
  }

  return cell;
}

}  // end namespace spin
}  // end namespace axom

#endif  // AXOM_SPIN_UNIFORMGRID_HPP_
