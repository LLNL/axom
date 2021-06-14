// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file IndirectionPolicies.hpp
 *
 * \brief Defines several indirection policies for slam
 *
 * Indirection policies encompass the underlying storage for indirection buffers
 * for a SLAM set, relation or map. A valid indirection policy must support the
 * following interface:
 *   * [required]
 *   * type alias IndirectionResult -- the type of the result of an indirection
 *      (const/nonconst and ref/nonref)
 *   * indirection() : IntType  -- returns the value of the element after
 *     indirection
 *   * hasIndirection(): bool -- returns whether there is an indirection
 *     buffer
 *   * isValid() : bool -- indicates whether the Indirection policy of the set
 *      is valid
 *   * [optional]
 *     * operator(): IntType -- alternate accessor for indirection
 *     * data() : ElementType* -- allows direct access to the underlying buffer
 *       (when this exists)
 *
 * \note Slam's Sets, Relations and Maps are not responsible for
 *  allocating/deallocating their own memory
 */

#ifndef SLAM_POLICIES_INDIRECTION_H_
#define SLAM_POLICIES_INDIRECTION_H_

#include "axom/core/Macros.hpp"
#include "axom/slic/interface/slic.hpp"

namespace axom
{
namespace slam
{
namespace policies
{
/**
 * \name OrderedSet_Indirection_Policies
 * \brief A few default policies for the indirection of an OrderedSet
 */

/// \{

/**
 * \brief A policy class for sets with no indirection
 */
template <typename PositionType, typename ElementType>
struct NoIndirection
{
  using IndirectionResult = ElementType;
  using ConstIndirectionResult = const ElementType;
  using IndirectionBufferType = struct
  { };

  NoIndirection() { }

  // This empty .ctor exists to satisfy IndirectionPolicy interface
  NoIndirection(IndirectionBufferType*) { }

  inline IndirectionResult indirection(PositionType pos) const
  {
    return static_cast<ElementType>(pos);
  }

  inline IndirectionResult operator()(PositionType pos) const
  {
    return indirection(pos);
  }

  IndirectionBufferType* data() { return nullptr; }

  bool hasIndirection() const { return false; }
  inline bool isValid(PositionType, PositionType, PositionType, bool) const
  {
    return true;
  }
};

/**
 * \brief A policy class for sets with array-based indirection
 */
template <typename PositionType, typename ElementType>
struct ArrayIndirection
{
  using IndirectionResult = ElementType&;
  using ConstIndirectionResult = const ElementType&;

  using IndirectionBufferType = ElementType;

  ArrayIndirection(IndirectionBufferType* buf = nullptr) : m_arrBuf(buf) { }

  IndirectionBufferType*& data() { return m_arrBuf; }

  inline ConstIndirectionResult indirection(PositionType pos) const
  {
    SLIC_ASSERT_MSG(hasIndirection(),
                    "SLAM::Set:ArrayIndirection -- Tried to dereference "
                      << " a null array in an array based indirection set.");
    return m_arrBuf[pos];
  }

  inline IndirectionResult indirection(PositionType pos)
  {
    SLIC_ASSERT_MSG(hasIndirection(),
                    "SLAM::Set:ArrayIndirection -- Tried to dereference "
                      << " a null array in an array based indirection set.");
    return m_arrBuf[pos];
  }

  inline ConstIndirectionResult operator()(PositionType pos) const
  {
    return indirection(pos);
  }

  inline IndirectionResult operator()(PositionType pos)
  {
    return indirection(pos);
  }

  bool hasIndirection() const { return m_arrBuf != nullptr; }

  inline bool isValid(PositionType size,
                      PositionType offset,
                      PositionType stride,
                      bool verboseOutput = false) const;

private:
  IndirectionBufferType* m_arrBuf;
};

/**
 * \brief A policy class for sets with stl vector-based indirection
 */
template <typename PositionType, typename ElementType>
struct STLVectorIndirection
{
  using IndirectionResult = ElementType&;
  using ConstIndirectionResult = const ElementType&;

  using VectorType = std::vector<ElementType>;
  using IndirectionBufferType = VectorType;

  STLVectorIndirection(IndirectionBufferType* buf = nullptr) : m_vecBuf(buf) { }

  IndirectionBufferType*& data() { return m_vecBuf; }
  IndirectionBufferType* const& data() const { return m_vecBuf; }

  inline ConstIndirectionResult indirection(PositionType pos) const
  {
    SLIC_ASSERT_MSG(hasIndirection(),
                    "SLAM::Set:STLVectorIndirection -- Tried to dereference "
                      << "a null vector in a vector based indirection set.");
    //SLIC_ASSERT_MSG( pos < m_vecBuf->size(),
    //  "SLAM::Set:STLVectorIndirection -- "
    //  << "Tried to access an out of bounds element at position "
    //  << pos << " in vector with only " << m_vecBuf->size() << " elements.");

    return (*m_vecBuf)[pos];
  }

  inline IndirectionResult indirection(PositionType pos)
  {
    SLIC_ASSERT_MSG(hasIndirection(),
                    "SLAM::Set:STLVectorIndirection -- Tried to dereference "
                      << "a null vector in a vector based indirection set.");

    return (*m_vecBuf)[pos];
  }

  inline IndirectionResult operator()(PositionType pos)
  {
    return indirection(pos);
  }

  inline ConstIndirectionResult operator()(PositionType pos) const
  {
    return indirection(pos);
  }

  bool hasIndirection() const { return m_vecBuf != nullptr; }

  inline bool isValid(PositionType size,
                      PositionType offset,
                      PositionType stride,
                      bool verboseOutput = false) const;

private:
  IndirectionBufferType* m_vecBuf;
};

/// \}

template <typename PosType, typename ElemType>
bool ArrayIndirection<PosType, ElemType>::isValid(PosType size,
                                                  PosType offset,
                                                  PosType stride,
                                                  bool verboseOutput) const
{
  AXOM_UNUSED_VAR(verboseOutput);

  // set of zero size is always valid
  if(size == 0) return true;

  bool bValid = true;

  // Check whether the set has elements, but the array ptr is null
  if(!hasIndirection())
  {
    SLIC_DEBUG_IF(verboseOutput,
                  "Array-based indirection set with non-zero size"
                    << " (size=" << size << ") requires a valid data buffer,"
                    << " but buffer pointer was null.");

    bValid = false;
  }
  else
  {
    // Check that none of the elements have negative indices within the array
    // Note: We do not have sufficient information about the array to know its
    // upper bound

    PosType firstEltInd = offset;
    PosType lastEltInd = (size - 1) * stride + offset;

    bool isRangeValid = (firstEltInd >= 0) && (lastEltInd >= 0);
    if(!isRangeValid)
    {
      SLIC_DEBUG_IF(
        verboseOutput,
        "Array-based indirection does not allow access "
          << "to data with lower addresses than its underlying pointer."
          << " Offset of " << offset << " leads to a first index of "
          << firstEltInd << "."
          << " Stride of " << stride << " and size of " << size
          << " leads to a last index of " << lastEltInd << ".");

      bValid = false;
    }
  }

  return bValid;
}

template <typename PosType, typename ElemType>
bool STLVectorIndirection<PosType, ElemType>::isValid(PosType size,
                                                      PosType offset,
                                                      PosType stride,
                                                      bool verboseOutput) const
{
  AXOM_UNUSED_VAR(verboseOutput);

  // always valid if set has zero size, even if indirection buffer is null
  if(size == 0) return true;

  bool bValid = true;

  // Otherwise, check whether the set has elements, but the array ptr is null
  if(!hasIndirection())
  {
    SLIC_DEBUG_IF(verboseOutput,
                  "Vector-based indirection set with non-zero size "
                    << "(size=" << size << ") requires a valid data buffer,"
                    << "but buffer pointer was null.");

    bValid = false;
  }
  else
  {
    // Verify underlying vector has sufficient storage for all set elements
    // Note: it is valid for the data buffer to have extra space
    PosType firstEltInd = offset;
    PosType lastEltInd = (size - 1) * stride + offset;
    PosType vecSize = static_cast<PosType>(m_vecBuf->size());

    bool isRangeValid = (0 <= firstEltInd) && (firstEltInd < vecSize) &&
      (0 <= lastEltInd) && (lastEltInd < vecSize);

    if(!isRangeValid)
    {
      SLIC_DEBUG_IF(
        verboseOutput,
        "Invalid vector-based IndirectionSet -- Data buffer "
          << "must be large enough to hold all elements of the set. "
          << "Underlying buffer size is " << vecSize << "."
          << " Offset of " << offset << " leads to a first index of "
          << firstEltInd << "."
          << " Stride of " << stride << " and size of " << size
          << " leads to a last index of " << lastEltInd << ".");

      bValid = false;
    }
  }

  return bValid;
}

}  // end namespace policies
}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_POLICIES_INDIRECTION_H_
