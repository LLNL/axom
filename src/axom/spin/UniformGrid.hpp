// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_UNIFORMGRID_HPP_
#define AXOM_SPIN_UNIFORMGRID_HPP_

#include "axom/core/utilities/Utilities.hpp"
#include "axom/core/execution/for_all.hpp"
#include "axom/core/Array.hpp"

#include "axom/slic/interface/slic.hpp"

#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/NumericArray.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/spin/RectangularLattice.hpp"

#include "axom/spin/policy/UniformGridStorage.hpp"

#if defined(AXOM_USE_RAJA)
  // RAJA includes
  #include "RAJA/RAJA.hpp"
#endif

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
 * \tparam ExecSpace The execution space in which operations will take place
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
template <typename T,
          int NDIMS,
          typename ExecSpace = axom::SEQ_EXEC,
          typename StoragePolicy = policy::DynamicGridStorage<T>>
class UniformGrid : StoragePolicy
{
public:
  static_assert((NDIMS == 3) || (NDIMS == 2),
                "Uniform grid dimensions must be 2 or 3.");

  /*! \brief The type used for specifying spatial extent of the contents */
  using BoxType = primal::BoundingBox<double, NDIMS>;

  /*! \brief The type used to query the index */
  using PointType = primal::Point<double, NDIMS>;

  struct QueryObject;

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
  UniformGrid(const double* lower_bound,
              const double* upper_bound,
              const int* res,
              int allocatorID = axom::execution_space<ExecSpace>::allocatorID());

  /*!
   * \brief Constructor specifying bounding box and number of bins.
   */
  UniformGrid(const BoxType& bbox,
              const int* res,
              int allocatorID = axom::execution_space<ExecSpace>::allocatorID());

  /*!
   * \brief Constructor specifying objects to initialize the UniformGrid with.
   */
  UniformGrid(const primal::NumericArray<int, NDIMS>& res,
              axom::ArrayView<const BoxType> bboxes,
              axom::ArrayView<const T> objs,
              int allocatorID = axom::execution_space<ExecSpace>::allocatorID());

  /*!
   * \brief Reinitializes a UniformGrid with an array of objects and associated
   *  bounding boxes.
   */
  void initialize(axom::ArrayView<const BoxType> bboxes,
                  axom::ArrayView<const T> objs);

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
  using StoragePolicy::getNumBins;

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
  using StoragePolicy::getBinContents;

  /*!
   * \brief Returns true if index is valid; that is, refers to a valid bin.
   *
   * The region of space indexed by this UniformGrid is divided into bins.
   * Each bin has a unique integer index.  If an integer is greater than or
   * equal to zero and less than or equal to the greatest bin index, this
   * method will return true.  Otherwise, the integer does not refer to any
   * bin and this method will return false.
   */
  using StoragePolicy::isValidIndex;

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
  using StoragePolicy::clear;

  /*!
   * \brief Inserts obj into each bin overlapped by BB.
   *
   * No error is signalled if BB falls partly or wholly outside the UniformGrid.
   *
   * \param [in] BB The region in which to record obj
   * \param [in] obj The object to insert into any bins overlapped by BB
   */
  void insert(const BoxType& BB, const T& obj);

  QueryObject getQueryObject() const;

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
   */
  static AXOM_HOST_DEVICE GridCell
  getClampedGridCell(const LatticeType& lattice,
                     const primal::NumericArray<int, NDIMS>& resolution,
                     const PointType& pt);

  /*! \brief Adds an object obj to the bin at index index */
  void addObj(const T& obj, int index);

  /*!
   * \brief Common constructor code for constructing rectangular lattice and
   *  bins in the uniform grid. Called after bounding box, resolution, and
   *  strides are set.
   */
  void initialize_grid();

private:
  BoxType m_boundingBox;
  LatticeType m_lattice;

  primal::NumericArray<int, NDIMS> m_resolution;
  primal::NumericArray<int, NDIMS> m_strides;

  DISABLE_COPY_AND_ASSIGNMENT(UniformGrid);
  DISABLE_MOVE_AND_ASSIGNMENT(UniformGrid);

};  //end class

template <typename T, int NDIMS, typename ExecSpace, typename StoragePolicy>
struct UniformGrid<T, NDIMS, ExecSpace, StoragePolicy>::QueryObject
  : StoragePolicy::ConstViewType
{
public:
  /*! \brief The type used for specifying spatial extent of the contents */
  using BoxType = primal::BoundingBox<double, NDIMS>;

  /*! \brief The type used to query the index */
  using PointType = primal::Point<double, NDIMS>;

  using ConstBinType =
    typename UniformGrid<T, NDIMS, ExecSpace, StoragePolicy>::ConstBinType;

  using LatticeType = RectangularLattice<NDIMS, double, int>;
  using GridCell = typename LatticeType::GridCell;

  QueryObject(BoxType bbox,
              LatticeType lattice,
              const primal::NumericArray<int, NDIMS>& resolution,
              const primal::NumericArray<int, NDIMS>& strides,
              const UniformGrid<T, NDIMS, ExecSpace, StoragePolicy>& from)
    : StoragePolicy::ConstViewType(from)
    , m_boundingBox(bbox)
    , m_lattice(lattice)
    , m_resolution(resolution)
    , m_strides(strides)
  { }

  /*!
   * \brief Returns the contents of the bin indicated by index.
   *
   * It is an error if index is invalid.
   * \param [in] index The index of the bin to retrieve.
   */
  using StoragePolicy::ConstViewType::getBinContents;

  AXOM_HOST_DEVICE ConstBinType getCandidates(const PointType& pt) const
  {
    if(!m_boundingBox.contains(pt))
    {
      return {};
    }

    // Find pt's integer grid cell within the uniform grid
    GridCell cell = UniformGrid::getClampedGridCell(m_lattice, m_resolution, pt);

    // compute the linear index (row-major) of the cell
    IndexType res = cell[0];
    for(int i = 1; i < NDIMS; ++i)
    {
      res += m_strides[i] * cell[i];
    }

    return getBinContents(res);
  }

  AXOM_HOST_DEVICE IndexType countCandidates(const BoxType& bbox) const
  {
    IndexType sumOfBinSizes = 0;
    auto lambdaCount = [&sumOfBinSizes, this](IndexType ibin) {
      sumOfBinSizes += getBinContents(ibin);
    };
    loopOverBboxes(bbox, lambdaCount);
    return sumOfBinSizes;
  }

  template <typename Func>
  AXOM_HOST_DEVICE void visitCandidates(const BoxType& bbox, Func&& evalFn) const
  {
    auto lambdaCandidates = [evalFn, this](IndexType ibin) {
      auto bin = getBinContents(ibin);
      for(IndexType ielem = 0; ielem < bin.size(); ielem++)
      {
        evalFn(bin[ielem]);
      }
    };
    loopOverBboxes(bbox, lambdaCandidates);
  }

private:
  template <typename Func>
  AXOM_HOST_DEVICE void loopOverBboxes(const BoxType& bbox, Func&& func) const
  {
    if(!m_boundingBox.intersectsWith(bbox))
    {
      return;
    }

    const GridCell lowerCell =
      getClampedGridCell(m_lattice, m_resolution, bbox.getMin());
    const GridCell upperCell =
      getClampedGridCell(m_lattice, m_resolution, bbox.getMax());

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
          func(i + jOffset);
        }
      }
    }
  }

  BoxType m_boundingBox;
  LatticeType m_lattice;

  primal::NumericArray<int, NDIMS> m_resolution;
  primal::NumericArray<int, NDIMS> m_strides;
};

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
template <typename T, int NDIMS, typename ExecSpace, typename StoragePolicy>
UniformGrid<T, NDIMS, ExecSpace, StoragePolicy>::UniformGrid(
  const double* lower_bound,
  const double* upper_bound,
  const int* res,
  int allocatorID)
  : UniformGrid(BoxType {PointType {lower_bound}, PointType {upper_bound}},
                res,
                allocatorID)
{ }

//------------------------------------------------------------------------------
template <typename T, int NDIMS, typename ExecSpace, typename StoragePolicy>
UniformGrid<T, NDIMS, ExecSpace, StoragePolicy>::UniformGrid(const BoxType& bbox,
                                                             const int* res,
                                                             int allocatorID)
  : StoragePolicy(allocatorID)
  , m_boundingBox(bbox)
  , m_resolution(-1)
{
  SLIC_ASSERT(res != nullptr);

  // set up the grid resolution
  m_resolution = primal::NumericArray<int, NDIMS>(res);

  initialize_grid();
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS, typename ExecSpace, typename StoragePolicy>
UniformGrid<T, NDIMS, ExecSpace, StoragePolicy>::UniformGrid(
  const primal::NumericArray<int, NDIMS>& res,
  axom::ArrayView<const BoxType> bboxes,
  axom::ArrayView<const T> objs,
  int allocatorID)
  : StoragePolicy(allocatorID)
  , m_resolution(res)
{
  initialize(bboxes, objs);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS, typename ExecSpace, typename StoragePolicy>
void UniformGrid<T, NDIMS, ExecSpace, StoragePolicy>::initialize_grid()
{
  // set up the row-major strides
  m_strides[0] = 1;
  for(int i = 1; i < NDIMS; ++i)
  {
    m_strides[i] = m_strides[i - 1] * m_resolution[i - 1];
  }

  // initialize space for the bins
  const int numBins = m_strides[NDIMS - 1] * m_resolution[NDIMS - 1];
  StoragePolicy::setNumBins(numBins);

  // scale the bounding box by a little to account for boundaries
  const double EPS = 1e-12;
  m_boundingBox.scale(1. + EPS);

  // set up the bounding box and lattice for point conversions
  m_lattice = rectangular_lattice_from_bounding_box(m_boundingBox, m_resolution);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS, typename ExecSpace, typename StoragePolicy>
void UniformGrid<T, NDIMS, ExecSpace, StoragePolicy>::initialize(
  axom::ArrayView<const BoxType> bboxes,
  axom::ArrayView<const T> objs)
{
  SLIC_ASSERT(bboxes.size() == objs.size());
  // get the global bounding box of all the objects
  BoxType global_box;
#ifdef AXOM_USE_RAJA
  using reduce_pol = typename axom::execution_space<ExecSpace>::reduce_policy;
  double infinity = std::numeric_limits<double>::max();
  double neg_infinity = std::numeric_limits<double>::lowest();
  RAJA::ReduceMin<reduce_pol, double> min_coord[NDIMS];
  RAJA::ReduceMax<reduce_pol, double> max_coord[NDIMS];
  for(int dim = 0; dim < NDIMS; dim++)
  {
    min_coord[dim].reset(infinity);
    max_coord[dim].reset(neg_infinity);
  }
  axom::for_all<ExecSpace>(
    bboxes.size(),
    AXOM_LAMBDA(IndexType idx) {
      for(int dim = 0; dim < NDIMS; dim++)
      {
        min_coord[dim].min(bboxes[idx].getMin()[dim]);
        max_coord[dim].max(bboxes[idx].getMax()[dim]);
      }
    });
  primal::Point<double, NDIMS> min_pt, max_pt;
  for(int dim = 0; dim < NDIMS; dim++)
  {
    min_pt[dim] = min_coord[dim].get();
    max_pt[dim] = max_coord[dim].get();
  }
  global_box = BoxType {min_pt, max_pt};
#else
  axom::for_all<ExecSpace>(bboxes.size(), [=, &global_box](IndexType idx) {
    global_box.addBox(bboxes[idx]);
  });
  m_boundingBox = global_box;
#endif

  // Now that we have the bounding box and resolution, initialize the
  // rectangular lattice and uniform grid storage.
  initialize_grid();

  const IndexType numBins = getNumBins();
  // 1. Get number of elements to insert into each bin
  axom::Array<IndexType> binCounts(numBins);
  // TODO: There's an error on operator[] if this isn't const and it only
  // happens for GCC 8.1.0
  const axom::ArrayView<IndexType> binCountsView = binCounts;

#ifdef AXOM_USE_RAJA
  using atomic_pol = typename axom::execution_space<ExecSpace>::atomic_policy;
#endif

  primal::NumericArray<int, NDIMS> strides = m_strides;
  axom::for_all<ExecSpace>(
    bboxes.size(),
    AXOM_LAMBDA(IndexType idx) {
      const BoxType bbox = bboxes[idx];
      const GridCell lowerCell =
        getClampedGridCell(m_lattice, m_resolution, bbox.getMin());
      const GridCell upperCell =
        getClampedGridCell(m_lattice, m_resolution, bbox.getMax());
      const int kLower = (NDIMS == 2) ? 0 : lowerCell[2];
      const int kUpper = (NDIMS == 2) ? 0 : upperCell[2];
      const int kStride = (NDIMS == 2) ? 1 : strides[2];

      for(IndexType k = kLower; k <= kUpper; ++k)
      {
        const IndexType kOffset = k * kStride;
        for(IndexType j = lowerCell[1]; j <= upperCell[1]; ++j)
        {
          const IndexType jOffset = j * strides[1] + kOffset;
          for(IndexType i = lowerCell[0]; i <= upperCell[0]; ++i)
          {
            const IndexType ibin = i + jOffset;
#ifdef AXOM_USE_RAJA
            RAJA::atomicAdd<atomic_pol>(&binCountsView[ibin], 1);
#else
            binCountsView[ibin]++;
#endif
          }
        }
      }
    });

  // 2. Resize bins with counts
  StoragePolicy::initialize(binCounts);

  // 3. Reset bin-specific counters
  binCounts.fill(0);

  // 4. Add elements to bins using a counting sort
  axom::for_all<ExecSpace>(
    bboxes.size(),
    AXOM_LAMBDA(IndexType idx) {
      const BoxType bbox = bboxes[idx];
      const GridCell lowerCell =
        getClampedGridCell(m_lattice, m_resolution, bbox.getMin());
      const GridCell upperCell =
        getClampedGridCell(m_lattice, m_resolution, bbox.getMax());
      const int kLower = (NDIMS == 2) ? 0 : lowerCell[2];
      const int kUpper = (NDIMS == 2) ? 0 : upperCell[2];
      const int kStride = (NDIMS == 2) ? 1 : strides[2];

      for(int k = kLower; k <= kUpper; ++k)
      {
        const int kOffset = k * kStride;
        for(int j = lowerCell[1]; j <= upperCell[1]; ++j)
        {
          const int jOffset = j * strides[1] + kOffset;
          for(int i = lowerCell[0]; i <= upperCell[0]; ++i)
          {
            const IndexType binIndex = i + jOffset;
            IndexType binCurrOffset = binCountsView[binIndex];
            StoragePolicy::get(binIndex, binCurrOffset) = objs[idx];
#ifdef AXOM_USE_RAJA
            RAJA::atomicAdd<atomic_pol>(&binCountsView[binIndex], 1);
#else
            binCountsView[binIndex]++;
#endif
          }
        }
      }
    });
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS, typename ExecSpace, typename StoragePolicy>
int UniformGrid<T, NDIMS, ExecSpace, StoragePolicy>::getBinIndex(
  const PointType& pt) const
{
  // Index is only valid when the point is within the bounding box
  if(!m_boundingBox.contains(pt))
  {
    return INVALID_BIN_INDEX;
  }

  // Find pt's integer grid cell within the uniform grid
  GridCell cell = getClampedGridCell(m_lattice, m_resolution, pt);

  // compute the linear index (row-major) of the cell
  int res = cell[0];
  for(int i = 1; i < NDIMS; ++i)
  {
    res += m_strides[i] * cell[i];
  }

  return res;
}

template <typename T, int NDIMS, typename ExecSpace, typename StoragePolicy>
bool UniformGrid<T, NDIMS, ExecSpace, StoragePolicy>::isBinEmpty(int index) const
{
  if(!isValidIndex(index))
  {
    return true;
  }

  return StoragePolicy::getBinContents(index).size() == 0;
}

//------------------------------------------------------------------------------

template <typename T, int NDIMS, typename ExecSpace, typename StoragePolicy>
const std::vector<int> UniformGrid<T, NDIMS, ExecSpace, StoragePolicy>::getBinsForBbox(
  const BoxType& BB) const
{
  std::vector<int> retval;

  if(!m_boundingBox.intersectsWith(BB))
  {
    return retval;
  }

  const GridCell lowerCell =
    getClampedGridCell(m_lattice, m_resolution, BB.getMin());
  const GridCell upperCell =
    getClampedGridCell(m_lattice, m_resolution, BB.getMax());

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
template <typename T, int NDIMS, typename ExecSpace, typename StoragePolicy>
void UniformGrid<T, NDIMS, ExecSpace, StoragePolicy>::addObj(const T& obj,
                                                             int index)
{
  SLIC_CHECK(isValidIndex(index));

  if(isValidIndex(index))
  {
    StoragePolicy::insert(index, obj);
  }
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS, typename ExecSpace, typename StoragePolicy>
void UniformGrid<T, NDIMS, ExecSpace, StoragePolicy>::insert(const BoxType& BB,
                                                             const T& obj)
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
template <typename T, int NDIMS, typename ExecSpace, typename StoragePolicy>
typename UniformGrid<T, NDIMS, ExecSpace, StoragePolicy>::QueryObject
UniformGrid<T, NDIMS, ExecSpace, StoragePolicy>::getQueryObject() const
{
  return QueryObject {m_boundingBox, m_lattice, m_resolution, m_strides, *this};
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS, typename ExecSpace, typename StoragePolicy>
typename UniformGrid<T, NDIMS, ExecSpace, StoragePolicy>::GridCell
UniformGrid<T, NDIMS, ExecSpace, StoragePolicy>::getClampedGridCell(
  const LatticeType& lattice,
  const primal::NumericArray<int, NDIMS>& resolution,
  const PointType& pt)
{
  GridCell cell = lattice.gridCell(pt);

  // clamp the grid coordinates to lie in a cell of the grid
  for(int i = 0; i < NDIMS; ++i)
  {
    cell[i] = axom::utilities::clampVal(cell[i], 0, resolution[i] - 1);
  }

  return cell;
}

}  // end namespace spin
}  // end namespace axom

#endif  // AXOM_SPIN_UNIFORMGRID_HPP_
