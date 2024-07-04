// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_STRIDED_STRUCTURED_INDEXING_HPP_
#define AXOM_MIR_STRIDED_STRUCTURED_INDEXING_HPP_

#include "axom/core/StackArray.hpp"
#include "axom/core/ArrayView.hpp"

namespace axom
{
namespace mir
{
namespace views
{

/**
 * \accelerated
 * \class StridedStructuredIndexing
 *
 * \brief This class encapsulates data for strided structured indexing and provides methods for creating/manipulating indices.
 *
 * \tparam NDIMS The number of dimensions.
 */
template <typename IndexT, int NDIMS = 3>
class StridedStructuredIndexing
{
public:
  using IndexType = IndexT;
  using LogicalIndex = axom::StackArray<axom::IndexType, NDIMS>;

  AXOM_HOST_DEVICE constexpr static int dimensions() { return NDIMS; }

  /**
   * \brief constructor
   */
  AXOM_HOST_DEVICE
  StridedStructuredIndexing() : m_dimensions(), m_offsets(), m_strides()
  {
    for(int i = 0; i < NDIMS; i++)
    {
      m_dimensions[i] = 1;
      m_offsets[i] = 0;
      m_strides[i] = 1;
    }
  }

  /**
   * \brief Constructor
   *
   * \param dims The number of zones in each logical dimension.
   * \param offsets The offset of the first zone from the mesh origin in each logical dimension.
   * \param strides The amount to stride when moving to the next element for each logical dimension.
   */
  AXOM_HOST_DEVICE
  StridedStructuredIndexing(const LogicalIndex &dims, const LogicalIndex &offsets, const LogicalIndex &strides) : m_dimensions(dims), m_offsets(offsets), m_strides(strides)
  {
  }

  /**
   * \brief Return the number of values in the index space.
   *
   * \return The number of values in the index space.
   */
  AXOM_HOST_DEVICE
  IndexType size() const
  {
     IndexType sz = 1;
     for(int i = 0; i < NDIMS; i++)
       sz *= m_dimensions[i];
     return sz;
  }

  /**
   * \brief Return the logical dimensions.
   *
   * \return The logical dimensions.
   */
  AXOM_HOST_DEVICE
  const LogicalIndex &logicalDimensions() const { return m_dimensions; }

  /**
   * \brief Return the j stride.
   *
   * \return The j stride to move up a row.
   */
  template <size_t _ndims = NDIMS>
  AXOM_HOST_DEVICE
  typename std::enable_if<_ndims >= 2, IndexType>::type
  jStride() const
  {
    return m_strides[1];
  }

  /**
   * \brief Return the k stride.
   *
   * \return The k stride to move forward a "page".
   */
  template <size_t _ndims = NDIMS>
  AXOM_HOST_DEVICE
  typename std::enable_if<_ndims == 3, IndexType>::type
  kStride() const
  {
    return m_strides[2];
  }

  /**
   * \brief Turn an index into a logical index.
   *
   * \param index The index to convert.
   *
   * \return The logical index that corresponds to the \a index.
   */
  /// @{

  template <size_t _ndims = NDIMS>
  AXOM_HOST_DEVICE
  typename std::enable_if<_ndims == 1, LogicalIndex>::type
  IndexToLogicalIndex(IndexType index) const
  {
    LogicalIndex logical;
    logical[0] = index - m_offsets[0];
    return logical;
  }

  template <size_t _ndims = NDIMS>
  AXOM_HOST_DEVICE
  typename std::enable_if<_ndims == 2, LogicalIndex>::type
  IndexToLogicalIndex(IndexType index) const
  {
    LogicalIndex logical;
    logical[0] = index % m_strides[1] - m_offsets[0];
    logical[1] = index / m_strides[1] - m_offsets[1];
    return logical;
  }

  template <size_t _ndims = NDIMS>
  AXOM_HOST_DEVICE
  typename std::enable_if<_ndims == 3, LogicalIndex>::type
  IndexToLogicalIndex(IndexType index) const
  {
    LogicalIndex logical;
    logical[0] = (index % m_strides[1]) - m_offsets[0];
    logical[1] = ((index % m_strides[2]) / m_strides[1]) - m_offsets[1];
    logical[2] = (index / m_strides[2]) - m_offsets[2];
    return logical;
  }

  /// @}

  /**
   * \brief Turn a logical index into a flat index.
   *
   * \param logical The logical indexto convert to a flat index.
   *
   * \return The index that corresponds to the \a logical index.
   */
  /// @{
  template <size_t _ndims = NDIMS>
  AXOM_HOST_DEVICE
  typename std::enable_if<_ndims == 1, IndexType>::type
  LogicalIndexToIndex(const LogicalIndex &logical) const
  {
    return ((m_offsets[0] + logical[0]) * m_strides[0]);
  }

  template <size_t _ndims = NDIMS>
  AXOM_HOST_DEVICE
  typename std::enable_if<_ndims == 2, IndexType>::type
  LogicalIndexToIndex(const LogicalIndex &logical) const
  {
    return ((m_offsets[0] + logical[0]) * m_strides[0]) +
           ((m_offsets[1] + logical[1]) * m_strides[1]);
  }

  template <size_t _ndims = NDIMS>
  AXOM_HOST_DEVICE
  typename std::enable_if<_ndims == 3, IndexType>::type
  LogicalIndexToIndex(const LogicalIndex &logical) const
  {
    return ((m_offsets[0] + logical[0]) * m_strides[0]) +
           ((m_offsets[1] + logical[1]) * m_strides[1]) +
           ((m_offsets[2] + logical[2]) * m_strides[2]);
  }

  /// @}

  /**
   * \brief Determines whether the indexing contains the supplied logical index.
   *
   * \param logical The logical index being tested.
   *
   * \return True if the logical index is within the index, false otherwise.
   */
  AXOM_HOST_DEVICE
  bool contains(const LogicalIndex &logical) const
  {
    bool retval = true;
    for(int i = 0; i < dimensions(); i++)
    {
      retval &= (logical[i] >= 0 && logical[i] < m_dimensions[i]);
    }
    return retval;
  }

  /**
   * \brief Determines whether the indexing contains the supplied index.
   *
   * \param index The index being tested.
   *
   * \return True if the index is within the index, false otherwise.
   */
  AXOM_HOST_DEVICE
  bool contains(const IndexType index) const
  {
    return contains(IndexToLogicalIndex(index));
  }

  /**
   * \brief Expand the current StridedStructuredIndexing by one in each dimension.
   *
   * \return An expanded StridedStructuredIndexing.
   */
  AXOM_HOST_DEVICE
  StridedStructuredIndexing expand() const
  {
    StridedStructuredIndexing retval(*this);
    for(int i = 0; i < dimensions(); i++)
      retval.m_dimensions[i]++;
    
    return retval;
  }

  /**
   * \brief Expand the current StridedStructuredIndexing by one in each dimension.
   *
   * \return An expanded StridedStructuredIndexing.
   */
  /// @{
  template <size_t _ndims = NDIMS>
  AXOM_HOST_DEVICE
  typename std::enable_if<_ndims == 1, StridedStructuredIndexing>::type
  expand() const
  {
    StridedStructuredIndexing retval(*this);
    retval.m_dimensions[0]++;
    return retval;
  }

  template <size_t _ndims = NDIMS>
  AXOM_HOST_DEVICE
  typename std::enable_if<_ndims == 2, StridedStructuredIndexing>::type
  expand() const
  {
    StridedStructuredIndexing retval(*this);
    retval.m_dimensions[0]++;
    retval.m_dimensions[1]++;
    retval.m_strides[1]++;
    return retval;
  }

  template <size_t _ndims = NDIMS>
  AXOM_HOST_DEVICE
  typename std::enable_if<_ndims == 3, StridedStructuredIndexing>::type
  expand() const
  {
    StridedStructuredIndexing retval(*this);
    retval.m_dimensions[0]++;
    retval.m_dimensions[1]++;
    retval.m_dimensions[2]++;
    const auto nx = retval.m_strides[1];
    const auto ny = retval.m_strides[2] / nx;
    retval.m_strides[1] = nx + 1;
    retval.m_strides[2] = (ny + 1) * (nx + 1);
    return retval;
  }

  /// @}

private:
  LogicalIndex m_dimensions{};
  LogicalIndex m_offsets{};
  LogicalIndex m_strides{};
};

} // end namespace views
} // end namespace mir
} // end namespace axom

#endif
