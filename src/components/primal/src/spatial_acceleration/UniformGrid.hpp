/*
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

#ifndef UNIFORMGRID_HPP_
#define UNIFORMGRID_HPP_

#include "slic/slic.hpp"

#include "primal/BoundingBox.hpp"
#include "primal/Point.hpp"


// C/C++ includes
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>
#include <stdio.h>

namespace axom {
namespace primal {

/*!
 *******************************************************************************
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
 *******************************************************************************
 */
template< typename T, int NDIMS >
class UniformGrid
{
public:
  /*! \brief The type used for specifying spatial extent of the contents */
  typedef BoundingBox< double, NDIMS > BoxType;

  /*! \brief The type used to query the index */
  typedef Point< double, NDIMS > PointType;

public:

  /*!
   *****************************************************************************
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
   * The points in the set {(2., y), (x, 4.)} are not included in the
   * UniformGrid.
   *****************************************************************************
   */
  UniformGrid(const double * lower_bound, const double * upper_bound, const int * res);

  /*! \brief Destructor: present for symmetry with constructor */
  ~UniformGrid();

  /*!
   *****************************************************************************
   * \brief Returns the index of the bin containing the specified point.
   *
   * If the specified point is outside the grid, the function returns a
   * special value (INVALID_BIN_INDEX).  Note that the upper boundaries of a
   * UniformGrid are not in a bin.  Thus, if a UniformGrid divides the region
   * from (0., 0.) to (10., 10.) into 10 x 10 bins, the code
   *     getBinIndex(PointType::makePoint(10., 10.))
   * will return INVALID_BIN_INDEX,
   * \param [in] pt The point to query.
   *****************************************************************************
   */
  int getBinIndex(const PointType & pt) const;

  /*! \brief Returns the number of bins in this UniformGrid. */
  int getNumBins() const;

  /*!
   *****************************************************************************
   * \brief Returns whether the bin specified by index is empty.
   *
   * Returns true if index is invalid.
   * \param [in] index The index of the bin to test.
   *****************************************************************************
   */
  bool isBinEmpty(int index) const;

  /*!
   *****************************************************************************
   * \brief Returns the contents of the bin indicated by index.
   *
   * It is an error if index is invalid.
   * \param [in] index The index of the bin to retrieve.
   *****************************************************************************
   */
  std::vector<T>& getBinContents(int index);

  /*!
   *****************************************************************************
   * \brief Returns the contents of the bin indicated by index.
   *
   * It is an error if index is invalid.  This is the const version.
   * \param [in] index The index of the bin to retrieve.
   *****************************************************************************
   */
  const std::vector<T>& getBinContents(int index) const;

  /*!
   *****************************************************************************
   * \brief Clears the bin indicated by index.
   *
   * No-op if index is invalid.
   * \param [in] index The index of the bin to clear.
   *****************************************************************************
   */
  void clear(int index);

  /*!
   *****************************************************************************
   * \brief Inserts obj into each bin overlapped by BB.
   *
   * No error is signalled if BB falls partly or wholly outside the UniformGrid.
   *
   * \param [in] BB The region in which to record obj
   * \param [in] obj The object to insert into any bins overlapped by BB
   *****************************************************************************
   */
  void insert(const BoxType& BB, const T& obj);

  /*!
   *****************************************************************************
   * \brief A special value indicating any location not in the UniformGrid.
   *
   * Returned by getBinIndex() if that function is passed an exterior point.
   * An error is thrown if getBinContents() is passed INVALID_BIN_INDEX.
   * clear(INVALID_BIN_INDEX) is a no-op; isBinEmpty(INVALID_BIN_INDEX) returns
   * true.
   *****************************************************************************
   */
  enum {
    INVALID_BIN_INDEX = -1
  };

protected:

  /*! \brief The default constructor should not be used, so it is protected. */
  UniformGrid();
    
private:

  PointType m_origin;
  double m_spacing[NDIMS];
  int m_resolution[NDIMS];

  void addObj(const T& obj, int index);
  bool isValidIndex(int index) const;

  struct Bin {
    std::vector<T> ObjectArray;
    //other stuff
  };
  std::vector<Bin> m_bins;

  DISABLE_COPY_AND_ASSIGNMENT(UniformGrid);
  DISABLE_MOVE_AND_ASSIGNMENT(UniformGrid);

};//end class

} //end namespace primal
} //end namespace axom



//------------------------------------------------------------------------------
// UniformGrid implementation
//------------------------------------------------------------------------------
namespace axom {
namespace primal {

//------------------------------------------------------------------------------
template< typename T, int NDIMS >
UniformGrid< T, NDIMS >::UniformGrid()
{

  SLIC_ASSERT((NDIMS == 3) || (NDIMS == 2));
    
  size_t newsize = 1;
  for (int i=0; i<NDIMS; ++i) {
    m_origin[i] = 0;
    m_spacing[i] = 1.0;
    m_resolution[i] = 100;
    newsize *= 100;
  }

  m_bins.resize(newsize);
}

template< typename T, int NDIMS >
UniformGrid< T, NDIMS >::UniformGrid(const double * lower_bound,
                                     const double * upper_bound,
                                     const int * res)
{
  SLIC_ASSERT(lower_bound != AXOM_NULLPTR);
  SLIC_ASSERT(upper_bound != AXOM_NULLPTR);
  SLIC_ASSERT(res != AXOM_NULLPTR);
  SLIC_ASSERT((NDIMS == 3) || (NDIMS == 2));

  size_t newsize = 1;
  for (int i=0; i<NDIMS;++i) {
    m_origin[i] = lower_bound[i];

    SLIC_ASSERT(lower_bound[i] <= upper_bound[i]);
    SLIC_ASSERT(res[i] > 0 );
    m_spacing[i] = (upper_bound[i] - lower_bound[i]) / res[i];
    m_resolution[i] = res[i];
    newsize *= res[i];
  }

  m_bins.resize(newsize);
}

template< typename T, int NDIMS >
UniformGrid< T, NDIMS >::~UniformGrid()
{
}

//------------------------------------------------------------------------------
template< typename T, int NDIMS >
int UniformGrid<T, NDIMS>::getBinIndex(const PointType & pt) const
{
  SLIC_ASSERT((NDIMS == 3) || (NDIMS == 2));

  int retval = 0;
  for (int i = 0; i < NDIMS; ++i) {
    int tmp = static_cast<int>(floor((pt[i] - m_origin[i]) / m_spacing[i]));

    if (tmp < 0 || tmp >= m_resolution[i]) {
      return INVALID_BIN_INDEX;
    }

    int factor = 1;
    for (int j = 0; j < i; ++j) {
      factor *= m_resolution[j];
    }

    retval += tmp * factor;
  }

  return retval;
}

template< typename T, int NDIMS >
bool UniformGrid<T, NDIMS>::isValidIndex(int index) const
{
  return index >=0 && index < static_cast<int>(m_bins.size());
}

template< typename T, int NDIMS >
int UniformGrid<T, NDIMS>::getNumBins() const
{
  return m_bins.size();
}

template< typename T, int NDIMS >
bool UniformGrid<T, NDIMS>::isBinEmpty(int index) const
{
  if (!isValidIndex(index)) {
    return true;
  }

  return m_bins[index].ObjectArray.empty();
}

//------------------------------------------------------------------------------
template< typename T, int NDIMS >
std::vector<T>& UniformGrid<T, NDIMS>::getBinContents(int index)
{
  SLIC_ASSERT(isValidIndex(index));

  return m_bins[index].ObjectArray;
}

template< typename T, int NDIMS >
const std::vector<T>& UniformGrid<T, NDIMS>::getBinContents(int index) const
{
  SLIC_ASSERT(isValidIndex(index));

  return m_bins[index].ObjectArray;
}

//------------------------------------------------------------------------------
template< typename T, int NDIMS >
void UniformGrid<T, NDIMS>::clear(int index)
{
  if (isValidIndex(index)) {
    m_bins[index].ObjectArray.clear();
  }
}

//------------------------------------------------------------------------------
template< typename T, int NDIMS >
void UniformGrid<T, NDIMS>::addObj(const T& obj, int index)
{
  SLIC_CHECK(!isValidIndex(index));

  if (isValidIndex(index)) {
    m_bins[index].ObjectArray.push_back(obj);
  }
}

//------------------------------------------------------------------------------
template< typename T, int NDIMS>
void UniformGrid<T, NDIMS>::insert(const BoxType& BB,
                                  const T& obj)
{
  SLIC_ASSERT((NDIMS == 3) || (NDIMS == 2));
   
  PointType bmin = BB.getMin();
  PointType bmax = BB.getMax();

  // Clamp the input bounding-box to at most the UniformGrid bounding box
  for (int dim = 0; dim < NDIMS; ++dim) {
    if (bmin[dim] < m_origin[dim]) {
      bmin[dim] = m_origin[dim];
    }
    double dimlimit = m_origin[dim] + m_spacing[dim] * (m_resolution[dim] - 1);
    if (bmax[dim] > dimlimit) {
      bmax[dim] = dimlimit;
    }
  }

  int start = getBinIndex(bmin);
  int end = getBinIndex(bmax);

  // Initialize bincount.  Note that NDIMS == 2 or NDIMS == 3.
  int bincount[3] = {1, 1, 1};

  // Guard against BB not overlapping the grid at all
  if (isValidIndex(start) && isValidIndex(end)) {
    int res = 1;
    // Find how many bboxes in each dimension (at least one)
    for (int dim = 0; dim < NDIMS; ++dim) {
      PointType extent(bmin);
      extent[dim] = bmax[dim];
      int bextent = getBinIndex(extent);
      bincount[dim] = 1 + (bextent - start) / res;

      // res accumulates the factor separating each dimension
      res *= m_resolution[dim];
    }

    const int x_res = m_resolution[0];

    const int y_res = m_resolution[1];
    const int kstep = x_res * y_res;

    for (int k = 0; k < bincount[2]; ++k) {
      const int k_offset = k * kstep;
      for (int j = 0; j < bincount[1]; ++j) {
        const int j_offset = j * x_res;
        for (int i = 0; i < bincount[0]; ++i) {
          addObj(obj, start + i + j_offset + k_offset);
        }
      }
    }
  }
}

}  /* end namespace primal */
}  /* end namespace axom */


#endif /* UNIFORMGRID_HPP_ */


