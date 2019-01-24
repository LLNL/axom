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

/**
 * \file IndirectionPolicies.hpp
 *
 * \brief Indirection policies for SLAM
 *
 * Indirection policies encompass the underlying storage for indirection buffers
 * for a SLAM set, relation or map. A valid indirection policy must support the
 * following interface:
 *   * [required]
 *   * typedef IndirectionResult -- the type of the result of an indirection
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
template<typename PositionType, typename ElementType>
struct NoIndirection
{
  typedef const ElementType IndirectionResult;
  typedef struct {}         IndirectionBufferType;

  NoIndirection() {}

  // This empty .ctor exists to satisfy IndirectionPolicy interface
  NoIndirection(IndirectionBufferType*) {}

  inline IndirectionResult          indirection(PositionType pos) const
  {
    return static_cast<ElementType>(pos);
  }

  inline IndirectionResult operator ()(PositionType pos)  const
  {
    return indirection(pos);
  }

  IndirectionBufferType* data() { return nullptr; }

  bool hasIndirection() const { return false; }
  inline bool isValid(PositionType, PositionType, PositionType, bool ) const
  {
    return true;
  }
};

/**
 * \brief A policy class for sets with array-based indirection
 */
template<typename PositionType, typename ElementType>
struct ArrayIndirection
{
  typedef const ElementType&  IndirectionResult;
  typedef ElementType IndirectionBufferType;

  ArrayIndirection(IndirectionBufferType* buf = nullptr)
    : m_arrBuf(buf) {}

  IndirectionBufferType*&   data() { return m_arrBuf; }

  inline IndirectionResult  indirection(PositionType pos) const
  {
    SLIC_ASSERT_MSG( hasIndirection(),
                     "SLAM::Set:ArrayIndirection -- Tried to dereference "
                     << " a null array in an array based indirection set.");
    return m_arrBuf[pos];
  }

  inline IndirectionResult operator ()(PositionType pos) const
  {
    return indirection(pos);
  }

  bool                              hasIndirection() const
  {
    return m_arrBuf != nullptr;
  }

  inline bool                       isValid(
    PositionType size,
    PositionType offset,
    PositionType stride,
    bool verboseOutput = false) const
  {
    // set of zero size is always valid
    if(size == 0)
      return true;

    bool bValid = true;

    // Check whether the set has elements, but the array ptr is null
    if( !hasIndirection() )
    {
      if(verboseOutput)
      {
        SLIC_DEBUG(
          "Array-based indirection set with non-zero size"
          << " (size=" << size << ") requires a valid data buffer,"
          << " but buffer pointer was null.");
      }

      bValid = false;
    }
    else
    {
      // Check that none of the elements have negative indices within the array
      // Note: We do not have sufficient information about the array to know its
      // upper bound

      PositionType firstEltInd = offset;
      PositionType lastEltInd = (size - 1) * stride + offset;

      bool isRangeValid = (firstEltInd >= 0) && (lastEltInd >= 0);
      if(!isRangeValid)
      {
        if(verboseOutput)
        {
          SLIC_DEBUG(
            "Array-based indirection does not allow access "
            << "to data with lower addresses than its underlying pointer."
            << " Offset of " << offset
            << " leads to a first index of " << firstEltInd << "."
            << " Stride of " << stride << " and size of " << size
            << " leads to a last index of " << lastEltInd << ".");
        }
        bValid = false;
      }
    }

    return bValid;
  }

private:
  IndirectionBufferType* m_arrBuf;
};

/**
 * \brief A policy class for sets with array-based indirection
 */
template<typename PositionType, typename ElementType>
struct STLVectorIndirection
{
  typedef std::vector<ElementType>  VectorType;
  typedef const ElementType&        IndirectionResult;
  typedef const VectorType IndirectionBufferType;


  STLVectorIndirection(IndirectionBufferType* buf = nullptr)
    : m_vecBuf(buf) {}

  IndirectionBufferType* &        data()       { return m_vecBuf; }
  IndirectionBufferType* const &  data() const { return m_vecBuf; }

  inline IndirectionResult        indirection(PositionType pos) const
  {
    SLIC_ASSERT_MSG(
      hasIndirection(),
      "SLAM::Set:STLVectorIndirection -- Tried to dereference "
      << "a null vector in a vector based indirection set.");
    //SLIC_ASSERT_MSG( pos < m_vecBuf->size(),
    //  "SLAM::Set:STLVectorIndirection -- "
    //  << "Tried to access an out of bounds element at position "
    //  << pos << " in vector with only " << m_vecBuf->size() << " elements.");

    return (*m_vecBuf)[pos];
  }
  inline IndirectionResult operator ()(PositionType pos) const
  {
    return indirection(pos);
  }

  bool                              hasIndirection() const
  {
    return m_vecBuf != nullptr;
  }

  inline bool                       isValid(
    PositionType size,
    PositionType offset,
    PositionType stride,
    bool verboseOutput = false) const
  {
    // If set has zero size, we are always valid (even if indirection buffer is
    // null)
    if(size == 0)
      return true;

    bool bValid = true;

    // Otherwise, check whether the set has elements, but the array ptr is null
    if( !hasIndirection() )
    {
      if(verboseOutput)
      {
        SLIC_DEBUG(
          "Vector-based indirection set with non-zero size (size="<< size <<
          ") requires a valid data buffer, but buffer pointer was null.");
      }

      bValid = false;
    }
    else
    {
      // Finally, check that the underlying vector has sufficient storage for
      // all set elements
      // Note that it is valid for the data buffer to have more space than the
      // set's positions
      PositionType firstEltInd = offset;
      PositionType lastEltInd = (size - 1) * stride + offset;
      PositionType vecSize = m_vecBuf->size();

      bool isRangeValid =
        (0 <= firstEltInd) && (firstEltInd < vecSize)
        && (0 <= lastEltInd) && (lastEltInd < vecSize);
      if(!isRangeValid)
      {
        if(verboseOutput)
        {
          SLIC_DEBUG(
            "Invalid vector-based IndirectionSet -- Data buffer "
            << "must be large enough to hold all elements of the set. "
            << "Underlying buffer size is " << vecSize << "."
            << " Offset of " << offset
            << " leads to a first index of " << firstEltInd << "."
            << " Stride of " << stride << " and size of " << size
            << " leads to a last index of " << lastEltInd << ".");
        }
        bValid = false;
      }
    }

    return bValid;
  }

private:
  IndirectionBufferType* m_vecBuf;
};

/// \}

} // end namespace policies
} // end namespace slam
} // end namespace axom

#endif // SLAM_POLICIES_INDIRECTION_H_
