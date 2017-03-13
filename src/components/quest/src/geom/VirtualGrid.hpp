/*
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

#ifndef VIRTUALGRID_HPP_
#define VIRTUALGRID_HPP_

#include "common/CommonTypes.hpp"
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


namespace quest {

/*!
 *******************************************************************************
 * \class VirtualGrid
 *
 * \brief A spatial index defined by origin, spacing, and resolution.
 *
 * \tparam T The type of object that this VirtualGrid will hold
 * \tparam NDIMS The VirtualGrid dimensionality: supported values are 2 or 3
 *
 * The VirtualGrid class is parameterized on the type of object it will contain
 * and the dimensionality of the grid.  Currently the VirtualGrid will index
 * 2D and 3D space.
 *
 * The VirtualGrid divides space into a number of bins, extending in a block
 * from an origin point and arranged in row-major order.  Each bin extends over
 * an axis-aligned rectangle or parallelepiped.  The size of a bin is specified
 * by the spacing constructor argument, and the number of bins is specified by
 * the res constructor argument.
 *
 * To use this class, a user will instantiate a VirtualGrid and insert a number
 * of objects (each with its own bounding box).  The insert operation puts the
 * object into some of the bins.  The user may retrieve the bin index of any
 * point, then retrieve any objects associated with the bin at that index.
 *******************************************************************************
 */
template< typename T, int NDIMS >
class VirtualGrid
{
public:
  /*! \brief The type used for specifying spatial extent of the contents */
  typedef axom::primal::BoundingBox< double, NDIMS > BoxType;

  /*! \brief The type used to query the index */
  typedef axom::primal::Point< double, NDIMS > PointType;

public:

  /*!
   *****************************************************************************
   * \brief Constructor specifying origin, size of bins, number of bins.
   *
   * The origin is specified as a point.  Each pointer argument is assumed to
   * point to an array of at least length NDIMS.
   *****************************************************************************
   */
  VirtualGrid(const PointType& origin, const double * spacing, const int * res);
  /*!
   *****************************************************************************
   * \brief Constructor specifying origin, size of bins, number of bins.
   *
   * The origin is specified as a double array.  Each pointer argument is
   * assumed to point to an array of at least length NDIMS.
   *****************************************************************************
   */
  VirtualGrid(const double * origin, const double * spacing, const int * res);

  /*! \brief Destructor: present for symmetry with constructor */
  ~VirtualGrid();

  /*!
   *****************************************************************************
   * \brief Returns the index of the bin containing the specified point.
   *
   * If the specified point is outside the grid, the function returns a
   * special value (INVALID_BIN_INDEX).
   * \param [in] pt The point to query.
   *****************************************************************************
   */
  int getBinIndex(const PointType & pt);

  /*! \brief Returns the number of bins in this VirtualGrid. */
  int getNumBins();

  /*!
   *****************************************************************************
   * \brief Returns whether the bin specified by index is empty.
   *
   * Returns true if index is invalid.
   * \param [in] index The index of the bin to test.
   *****************************************************************************
   */
  bool binEmpty(int index);

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
   * No error is signalled if BB falls partly or wholly outside the VirtualGrid.
   *
   * \param [in] BB The region in which to record obj
   * \param [in] obj The object to insert into any bins overlapped by BB
   *****************************************************************************
   */
  void insert(const BoxType& BB, const T& obj);

  /*!
   *****************************************************************************
   * \brief A special value indicating any location not in the VirtualGrid.
   *
   * Returned by getBinIndex() if that function is passed an exterior point.
   * An error is thrown if getBinContents() is passed INVALID_BIN_INDEX.
   * clear(INVALID_BIN_INDEX) is a no-op; binEmpty(INVALID_BIN_INDEX) returns
   * true.
   *****************************************************************************
   */
  enum {
    INVALID_BIN_INDEX = -1
  };

protected:

  /*! \brief The default constructor should not be used, so it is protected. */
  VirtualGrid();
    
private:

  PointType m_origin;
  double m_spacing[NDIMS];
  int m_resolution[NDIMS];

  void addObj(const T& obj, int index);
  bool isValidIndex(int index) const;
  void insertIntoBins(int start, int * binlength, const T& obj);

  struct Bin {
    std::vector<T> ObjectArray;
    //other stuff
  };
  std::vector<Bin> m_bins;

  DISABLE_COPY_AND_ASSIGNMENT(VirtualGrid);
  DISABLE_MOVE_AND_ASSIGNMENT(VirtualGrid);

};//end class

}//end namespace quest



//------------------------------------------------------------------------------
// VirtualGrid implementation
//------------------------------------------------------------------------------
namespace quest
{

//------------------------------------------------------------------------------
template< typename T, int NDIMS >
VirtualGrid< T, NDIMS >::VirtualGrid()
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
VirtualGrid< T, NDIMS >::VirtualGrid(const PointType& origin,
                                     const double * step,
                                     const int * res)
{
  SLIC_ASSERT(step != ATK_NULLPTR);
  SLIC_ASSERT(res != ATK_NULLPTR);
  SLIC_ASSERT((NDIMS == 3) || (NDIMS == 2));

  size_t newsize = 1;
  for (int i=0; i<NDIMS;++i) {
    m_origin[i] = origin[i];

    SLIC_ASSERT(step[i] !=0 );
    m_spacing[i] = step[i];
    m_resolution[i] = res[i];
    newsize *= res[i];
  }

  m_bins.resize(newsize);
}

template< typename T, int NDIMS >
VirtualGrid< T, NDIMS >::VirtualGrid(const double * origin,
                                     const double * step,
                                     const int * res)
{
  SLIC_ASSERT(origin != ATK_NULLPTR);
  SLIC_ASSERT(step != ATK_NULLPTR);
  SLIC_ASSERT(res != ATK_NULLPTR);
  SLIC_ASSERT((NDIMS == 3) || (NDIMS == 2));

  size_t newsize = 1;
  for (int i=0; i<NDIMS;++i) {
    m_origin[i] = origin[i];
    
    SLIC_ASSERT(step[i] !=0 );
    m_spacing[i] = step[i];
    m_resolution[i] = res[i];
    newsize *= res[i];
  }

  m_bins.resize(newsize);
}

template< typename T, int NDIMS >
VirtualGrid< T, NDIMS >::~VirtualGrid()
{
}

//------------------------------------------------------------------------------
template< typename T, int NDIMS >
int VirtualGrid<T, NDIMS>::getBinIndex(const PointType & pt)
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
bool VirtualGrid<T, NDIMS>::isValidIndex(int index) const
{
  return index >=0 && index < static_cast<int>(m_bins.size());
}

template< typename T, int NDIMS >
int VirtualGrid<T, NDIMS>::getNumBins()
{
  return m_bins.size();
}

template< typename T, int NDIMS >
bool VirtualGrid<T, NDIMS>::binEmpty(int index)
{
  if (!isValidIndex(index)) {
    return true;
  }

  return m_bins[index].ObjectArray.empty();
}

//------------------------------------------------------------------------------
template< typename T, int NDIMS >
std::vector<T>& VirtualGrid<T, NDIMS>::getBinContents(int index)
{
  SLIC_ASSERT(isValidIndex(index));

  return m_bins[index].ObjectArray;
}

template< typename T, int NDIMS >
const std::vector<T>& VirtualGrid<T, NDIMS>::getBinContents(int index) const
{
  SLIC_ASSERT(isValidIndex(index));

  return m_bins[index].ObjectArray;
}

//------------------------------------------------------------------------------
template< typename T, int NDIMS >
void VirtualGrid<T, NDIMS>::clear(int index)
{
  if (isValidIndex(index)) {
    m_bins[index].ObjectArray.clear();
  }
}

//------------------------------------------------------------------------------
template< typename T, int NDIMS >
void VirtualGrid<T, NDIMS>::addObj(const T& obj, int index)
{
  SLIC_CHECK(!isValidIndex(index));

  if (isValidIndex(index)) {
    m_bins[index].ObjectArray.push_back(obj);
  }
}


//------------------------------------------------------------------------------
template< typename T, int NDIMS>
void VirtualGrid<T, NDIMS>::insert(const BoxType& BB,
                                  const T& obj)
{
  SLIC_ASSERT((NDIMS == 3) || (NDIMS == 2));
   
  PointType bmin = BB.getMin();
  PointType bmax = BB.getMax();

  // Clamp the input bounding-box to at most the VirtualGrid bounding box
  for (int dim = 0; dim < NDIMS; ++dim) {
    if (bmin[dim] < m_origin[dim]) {
      bmin[dim] = m_origin[dim];
    }
    double dimlimit = m_origin[dim] + m_spacing[dim] * (m_resolution[dim] - 1);
    if (bmax[dim] > dimlimit) {
      bmax[dim] = dimlimit;
    }
  }

  int bincount[NDIMS];
  int start = getBinIndex(bmin);
  int end = getBinIndex(bmax);

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

#if NDIMS == 2
    for (int j = 0; j < bincount[1]; ++j) {
      const int j_offset = j * x_res;
      for (int i = 0; i < bincount[0]; ++i) {
        addObj(obj, start + i + j_offset);
      }
    }
#else
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
#endif
  }
}

}  /* end namespace quest */


#endif /* VIRTUALGRID_HPP_ */


