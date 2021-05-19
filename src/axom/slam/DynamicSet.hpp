// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file DynamicSet.hpp
 *
 * \brief Contains a DynamicSet class, whose size can change dynamically
 * at runtime
 */

#ifndef SLAM_DYNAMIC_SET_H_
#define SLAM_DYNAMIC_SET_H_

#include "axom/slam/OrderedSet.hpp"

namespace axom
{
namespace slam
{
/**
 * \class DynamicSet
 * \brief A Set class that supports dynamically adding and removing set items
 *
 * \detail An entry in the set is valid if it is not equal to INVALID_ENTRY.
 *
 * An example to traverse the elements
 * \code
 * DynamicSet<> some_set;
 * ... // initialize some_set
 *
 * const IndexType N = some_set.size()
 * for(IndexType i=0; i< N; ++i)
 * {
 *   if( some_set.isValidEntry(i) )
 *   {
 *     ElementType el = some_set[i];
 *
 *     ... // do something with el
 *
 *   } // END if the entry is valid
 * } //END for all set entries
 *
 * \endcode
 */

template <typename PosType = slam::DefaultPositionType,
          typename ElemType = slam::DefaultElementType,
          typename SizePolicy = policies::DynamicRuntimeSize<PosType>,
          typename OffsetPolicy = policies::ZeroOffset<PosType>,
          typename StridePolicy = policies::StrideOne<PosType>>
class DynamicSet : public Set<PosType, ElemType>,
                   SizePolicy,
                   OffsetPolicy,
                   StridePolicy
{
public:
  using PositionType = PosType;
  using ElementType = PosType;
  using SetVectorType = std::vector<ElementType>;
  using SizePolicyType = SizePolicy;

  enum
  {
    INVALID_ENTRY = ~0  ///< value to mark indices of deleted elements
  };

  struct SetBuilder;

public:
  /**
   * \brief Constructor for a DynamicSet
   * \param size The initial size of the set
   *
   * \note The set entries will be initialized such that set[i] = i
   */
  DynamicSet(PositionType size = SizePolicyType::DEFAULT_VALUE)
    : SizePolicy(size)
  {
    fill_array_default(size);
  };

  /**
   * \brief Constructor for a DynamicSet from a SetBuilder
   *
   * \note The set entries will be initialized such that set[i] = i
   */
  DynamicSet(const SetBuilder& builder)
    : SizePolicy(builder.m_size)
    , OffsetPolicy(builder.m_offset)
    , StridePolicy(builder.m_stride)
  {
    fill_array_default(builder.m_size.size());
  }

  //~DynamicSet();

public:
  /**
   * \class SetBuilder
   * \brief Helper class for constructing a DynamicSet.
   */
  struct SetBuilder
  {
    friend class DynamicSet;

    /** \brief Set the size of the DynamicSet using SizePolicy */
    SetBuilder& size(PositionType sz)
    {
      m_size = SizePolicy(sz);
      return *this;
    }

    /** \brief Set the offset of the DynamicSet using OffsetPolicy */
    SetBuilder& offset(PositionType off)
    {
      m_offset = OffsetPolicy(off);
      return *this;
    }

    /** \brief Set the stride of the DynamicSet using StridePolicy */
    SetBuilder& stride(PositionType str)
    {
      m_stride = StridePolicy(str);
      return *this;
    }

  private:
    SizePolicy m_size;
    OffsetPolicy m_offset;
    StridePolicy m_stride;
  };

public:
  /// \name DynamicSet element access functions
  /// @{

  /**
   * \brief Access the element at position \a pos
   *
   * \pre pos must be between 0 and size()
   */
  ElementType at(PositionType pos) const { return operator[](pos); };

  /**
   * \brief Access the element at position \a pos
   *
   * \pre pos must be between 0 and size()
   */
  ElementType operator[](IndexType pos) const
  {
    verifyPosition(pos);
    return m_data[pos];
  };

  /**
   * \brief Access the element at position \a pos
   *
   * \pre pos must be between 0 and size()
   */
  ElementType& operator[](IndexType pos)
  {
    verifyPosition(pos);
    return m_data[pos];
  };

  /** \brief Returns a reference to the underlying set data */
  SetVectorType& data() { return m_data; }

  /** \brief Returns a const reference to the underlying set data */
  const SetVectorType& data() const { return m_data; }

  /**
   * \brief Given a value, find the index of the first entry containing it
   *
   * \return The index of the first element with value \a e, or INVALID_ENTRY
   * if none can be found.
   * \note This is an O(n) operation
   */
  IndexType findIndex(ElementType e)
  {
    for(unsigned int i = 0; i < m_data.size(); ++i)
    {
      if(m_data[i] == e) return i;
    }
    return INVALID_ENTRY;
  };

  /// @}

public:
  /// \name Functions that deal with the set cardinality
  /// @{

  /**
   * \brief Returns the number of possible elements in the set
   *
   * \note Not all elements are necessarily valid since some elements
   * could have been deleted
   * \sa numberOfValidEntries(), isValidEntry()
   */
  PositionType size() const
  {
    return static_cast<PositionType>(m_data.size());
  };

  /** \brief Uses \a SizePolicy::empty() to determine if the set is empty */
  bool empty() const { return SizePolicy::empty(); };

  /**
   * \brief Return the number of valid entries in the set.
   *
   * \detail This is an O(n) operation, because the class makes no assumption
   * that data was not changed by the user
   */
  PositionType numberOfValidEntries() const
  {
    PositionType nvalid = 0;

    const int sz = static_cast<int>(m_data.size());
    for(int i = 0; i < sz; ++i)
    {
      nvalid += (m_data[i] != INVALID_ENTRY) ? 1 : 0;
    }
    return nvalid;
  }

  /** \brief Returns true if this set is a subset of another set */
  bool isSubset() const { return false; };

  /// @}

public:
  /// \name Functions that deal with validity checks
  /// @{

  /**
   * \brief Predicate to check if the entry at index \a i is valid
   *
   * The entry is valid when 0 <= i < size() and the value at index
   * \a i is not marked as \a INVALID_ENTRY
   */
  bool isValidEntry(IndexType i) const
  {
    return i >= 0 && i < static_cast<IndexType>(m_data.size()) &&
      m_data[i] != INVALID_ENTRY;
  };

  /**
   * \brief Returns true if the DynamicSet instance is valid
   *
   * A DynamicSet is valid if each of its policies claim it to be valid.
   * This includes its \a SizePolicy, \a OffsetPolicy and \a StridePolicy
   */
  bool isValid(bool verboseOutput = false) const
  {
    bool bValid = SizePolicy::isValid(verboseOutput) &&
      OffsetPolicy::isValid(verboseOutput) &&
      StridePolicy::isValid(verboseOutput);

    return bValid;
  };

  /// @}

public:
  /// \name Functions that modify the set cardinality
  /// @{

  /**
   * \brief Insert an entry at the end of the set with value = ( size()-1 )
   */
  IndexType insert() { return insert(static_cast<IndexType>(m_data.size())); }

  /**
   * \brief Insert an entry at the end of the set with the given value.
   * \param val the value of the inserted entry
   */
  IndexType insert(ElementType val)
  {
    m_data.push_back(val);
    return static_cast<IndexType>(m_data.size() - 1);
  };

  /**
   * \brief Mark the corresponding entry as invalid
   *
   * \note It is not a problem to mark an INVALID_ENTRY as INVALID_ENTRY
   */
  void remove(IndexType idx)
  {
    verifyPosition(idx);

    if(m_data[idx] != INVALID_ENTRY)
    {
      m_data[idx] = INVALID_ENTRY;
    }
  };

  /// @}

private:
  /** \brief Debug check that the index \a pos is not out-of-range */
  void verifyPosition(PositionType AXOM_DEBUG_PARAM(pos)) const
  {
    SLIC_ASSERT_MSG(
      (pos >= 0) && (pos < static_cast<PositionType>(m_data.size())),
      "SLAM::DynamicSet -- requested out-of-range element at position "
        << pos << ", but set only has " << m_data.size() << " elements.");
  };

  /** Fill each entry of the set such that its value is equal to its index. */
  void fill_array_default(PositionType size)
  {
    if(size < 0) return;

    m_data.resize(size);
    for(int i = 0; i < size; ++i)
    {
      m_data[i] = i;
    }
  }

private:
  SetVectorType m_data;
};

}  // end namespace slam
}  // end namespace axom

#endif
