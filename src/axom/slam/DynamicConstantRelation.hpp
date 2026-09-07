// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/**
 * \file DynamicConstantRelation.hpp
 *
 * \brief API for a topological relation between two sets in which entities
 * from the first set can be related to a constant number of entities from
 * the second set. For example, in a triangle mesh, each triangle is
 * incident to three vertices.
 *
 * Connectivity can change at runtime.
 */

#include "axom/config.hpp"
#include "axom/slic.hpp"

#include "axom/slam/Set.hpp"
#include "axom/slam/Relation.hpp"
#include "axom/slam/OrderedSet.hpp"
#include "axom/slam/DynamicSet.hpp"
#include "axom/slam/policies/CardinalityPolicies.hpp"
#include "axom/slam/policies/PolicyTraits.hpp"

#include "axom/fmt.hpp"

#include <vector>

namespace axom::slam
{
/**
 * \class DynamicConstantRelation
 * \brief  A relation class with constant cardinality that supports
 * adding, removing and modifying set relations.
 *
 * A DynamicConstantRelation encodes the relation between two sets,
 * A FromSet and a ToSet, where the cardinality of the relation from
 * each element of the FromSet to the ToSet is fixed to a constant value.
 * For example, each triangle in the triangle set of a triangle mesh
 * has three incident vertices from set of vertices.
 *
 * Connectivity stores to-set positions, not to-set element values. An unused
 * entry contains INVALID_INDEX. isValidEntry(from) checks that the from-set
 * entry is valid and at least one associated position differs from INVALID_INDEX.
 *
 * The relation owns its connectivity vector and borrows its sets. Those sets
 * must outlive the relation. Use updateSizes() after growing the from-set.
 */
template <typename PosType,   //= slam::DefaultPositionType,
          typename ElemType,  // = slam::DefaultElementType,
          typename CardinalityPolicy>
class DynamicConstantRelation : public /*Relation,*/ CardinalityPolicy
{
public:
  enum
  {
    INVALID_INDEX = ~0  ///< value to mark indices of deleted elements
  };

  using FromSetType = DynamicSet<PosType, ElemType>;
  using ToSetType = DynamicSet<PosType, ElemType>;
  using FromPositionType = typename FromSetType::PositionType;
  using ToPositionType = typename ToSetType::PositionType;

  using RelationVec = std::vector<ToPositionType>;

  using BeginsSizePolicy = typename CardinalityPolicy::RelationalOperatorSizeType;

  using STLIndirection = policies::STLVectorIndirection<FromPositionType, ToPositionType>;
  using RelationSubset = OrderedSet<FromPositionType,
                                    ToPositionType,
                                    BeginsSizePolicy,
                                    policies::RuntimeOffset<FromPositionType>,
                                    policies::StrideOne<FromPositionType>,
                                    STLIndirection>;

  // types for iterator
  using RelationIterator = typename RelationSubset::iterator;
  using RelationIteratorPair = typename RelationSubset::iterator_pair;

  using RelationConstIterator = typename RelationSubset::const_iterator;
  using RelationConstIteratorPair = typename RelationSubset::const_iterator_pair;

public:
  /// \brief Default constructor with empty set for toSet and fromSet
  DynamicConstantRelation()
    : m_fromSet(policies::EmptySetTraits<FromSetType>::emptySet())
    , m_toSet(policies::EmptySetTraits<ToSetType>::emptySet())
  { }

  /// \brief Construct a DynamicConstantRelation from the given \a fromSet to \a toSet
  DynamicConstantRelation(FromSetType* fromSet, ToSetType* toSet)
    : CardinalityPolicy(policies::EmptySetTraits<FromSetType>::isEmpty(fromSet) ? 0 : fromSet->size())
    , m_fromSet(fromSet)
    , m_toSet(toSet)
    , m_currentFromSize(fromSet == nullptr ? 0 : m_fromSet->size())
  {
    updateSizes();
  };

public:
  /// \name DynamicConstantRelation set accessors
  /// @{

  /// \brief Returns a pointer to the relation's from-set.
  FromSetType* fromSet() { return m_fromSet; }

  /// \overload
  const FromSetType* fromSet() const { return m_fromSet; }

  /// \brief Returns a pointer to the relation's to-set.
  ToSetType* toSet() { return m_toSet; }

  /// \overload
  const ToSetType* toSet() const { return m_toSet; }

  /// @}

public:
  /// \name DynamicConstantRelation iterator interface
  /// @{

  /**
   * \brief Begin iterating over the to-set positions related to fromSetInd.
   * \pre 0 <= fromSetInd < fromSet()->size()
   */
  RelationIterator begin(FromPositionType fromSetInd)
  {
    verifyPosition(fromSetInd);
    return (*this)[fromSetInd].begin();
  }

  /**
   * \brief Begin const iteration over the to-set positions related to fromSetInd.
   * \pre 0 <= fromSetInd < fromSet()->size()
   */
  RelationConstIterator begin(FromPositionType fromSetInd) const
  {
    verifyPosition(fromSetInd);
    return (*this)[fromSetInd].begin();
  }

  /**
   * \brief Return an iterator past the to-set positions related to fromSetInd.
   * \pre 0 <= fromSetInd < fromSet()->size()
   */
  RelationIterator end(FromPositionType fromSetInd)
  {
    verifyPosition(fromSetInd);
    return (*this)[fromSetInd].end();
  }

  /**
   * \brief Return a const iterator past the to-set positions related to fromSetInd.
   * \pre 0 <= fromSetInd < fromSet()->size()
   */
  RelationConstIterator end(FromPositionType fromSetInd) const
  {
    verifyPosition(fromSetInd);
    return (*this)[fromSetInd].end();
  }

  /**
   * \brief Return begin and end iterators over the to-set positions related to fromSetInd.
   * \pre 0 <= fromSetInd < fromSet()->size()
   */
  RelationIteratorPair range(FromPositionType fromSetInd) { return (*this)[fromSetInd].range(); }

  /**
   * \brief Return const begin and end iterators over the to-set positions related to fromSetInd.
   * \pre 0 <= fromSetInd < fromSet()->size()
   */
  RelationConstIteratorPair range(FromPositionType fromSetInd) const
  {
    return (*this)[fromSetInd].range();
  }

  /// @}

public:
  /// \name DynamicConstantRelation per-element relation access functions
  /// @{
  ///

  /**
   * \brief Return the to-set positions associated with fromSetIndex.
   * \param fromSetIndex The index of an element in the FromSet
   */
  RelationSubset const at(FromPositionType fromSetIndex) const
  {
    verifyPosition(fromSetIndex);
    return operator[](fromSetIndex);
  }

  RelationSubset at(FromPositionType fromSetIndex)
  {
    verifyPosition(fromSetIndex);
    return operator[](fromSetIndex);
  }

  /**
   * \brief Return the to-set positions associated with fromSetIndex.
   * \param fromSetIndex The index of an element in the FromSet
   * \note This function does not grow the relation. Use updateSizes() after
   *       growing the from-set. insert() and modify() can also expand storage.
   */
  RelationSubset const operator[](FromPositionType fromSetIndex) const
  {
    // NOTE: Need to const_cast the pointer to the vector
    // since SetBuilder, and the IndirectionPolicy don't
    // currently support const buffers
    // TODO: Fix this!

    verifyPosition(fromSetIndex);
    using SetBuilder = typename RelationSubset::SetBuilder;
    return SetBuilder()
      .size(relationCardinality())
      //.offset( CardinalityPolicy::offset( fromSetIndex) )
      .offset(fromSetIndex * relationCardinality())
      .data(const_cast<RelationVec*>(&m_relationsVec));
  }

  RelationSubset operator[](FromPositionType fromSetIndex)
  {
    verifyPosition(fromSetIndex);
    using SetBuilder = typename RelationSubset::SetBuilder;
    return SetBuilder()
      .size(relationCardinality())
      //.offset( CardinalityPolicy::offset( fromSetIndex) )
      .offset(fromSetIndex * relationCardinality())
      .data(&m_relationsVec);
  }

  /**
   * \brief Returns the cardinality of the set of entities in the ToSet
   * related to the element with index \a fromSetIndex in the FromSet
   * \param fromSetIndex The index of an element in the FromSet
   */
  FromPositionType size(FromPositionType fromSetIndex) const
  {
    verifyPosition(fromSetIndex);
    return relationCardinality();
  }

  /// @}

  /// \brief Returns the cardinality of the FromSet
  inline FromPositionType size() const { return m_currentFromSize; }

public:
  /// \name DynamicConstantRelation validity check functions
  /// @{
  ///

  /**
   * \brief Returns the number of valid entries in the FromSet
   *
   * An element of the FromSet is considered valid with respect to a
   * DynamicConstantRelation when it is valid in the FromSet and when
   * its relation set is not marked as invalid.
   * \sa isValidEntry()
   */
  FromPositionType numberOfValidEntries() const
  {
    FromPositionType nvalid = 0;
    const FromPositionType N = size();
    for(FromPositionType i = 0; i < N; ++i)
    {
      nvalid += isValidEntry(i);
    }
    return nvalid;
  }

  /**
   * \brief Check whether idx is valid in the from-set and has any assigned connectivity.
   * \details At least one associated value must differ from INVALID_INDEX.
   *          This does not check the bounds of those to-set positions.
   */
  bool isValidEntry(FromPositionType idx) const
  {
    if(m_fromSet->isValidEntry(idx))
    {
      const auto SZ = relationCardinality();
      const auto beg_idx = idx * SZ;
      for(auto i = beg_idx; i < (beg_idx + SZ); ++i)
      {
        if(m_relationsVec[i] != INVALID_INDEX)
        {
          return true;
        }
      }
    }
    return false;
  }

  /**
   * \brief Predicate to check if the DynamicConstantRelation instance is valid
   */
  bool isValid(bool verboseOutput = false) const;

  /// @}

public:
  /// \name DynamicConstantRelation functions that modify the relation
  /// @{
  ///

  /**
   * \brief Inserts a new entry into the relation at the first invalid index
   * \param fromSetIndex The index of the element in the FromSet
   * \param toSetIndex The index of the element in the ToSet
   * to associate with \a fromSetIndex
   */
  void insert(FromPositionType fromSetIndex, FromPositionType toSetIndex)
  {
    expandSizeIfNeeded(fromSetIndex + 1);
    verifyPosition(fromSetIndex);

    //find the first invalid place to put it
    const auto SZ = relationCardinality();
    for(int i = 0; i < SZ; ++i)
    {
      const auto idx = SZ * fromSetIndex + i;
      if(m_relationsVec[idx] == INVALID_INDEX)
      {
        m_relationsVec[idx] = toSetIndex;
        return;
      }
    }

    //The entry was not inserted
    SLIC_WARNING("Relation from " << fromSetIndex << " to " << toSetIndex
                                  << " was not inserted because the entry is full.");
  }

  /**
   * \brief Assign a to-set position at a local offset for fromSetIndex.
   *
   * \param fromSetIndex Position in the from-set.
   * \param offset Position within the related subset.
   * \param toSetIndex Position in the to-set.
   * \pre 0 <= offset < relationCardinality()
   * \note Expands storage if needed. For existing entries, the assignment is
   *       equivalent to relation[fromSetIndex][offset] = toSetIndex.
   */
  void modify(FromPositionType fromSetIndex, FromPositionType offset, FromPositionType toSetIndex)
  {
    expandSizeIfNeeded(fromSetIndex + 1);
    verifyPosition(fromSetIndex);
    m_relationsVec[relationCardinality() * fromSetIndex + offset] = toSetIndex;
  }

  /// \brief Mark all values in entry \a fromSetIndex as invalid.
  void remove(FromPositionType fromSetIndex)
  {
    if(!isValidEntry(fromSetIndex))
    {
      return;
    }

    const auto SZ = relationCardinality();
    for(int i = 0; i < SZ; ++i)
    {
      m_relationsVec[SZ * fromSetIndex + i] = INVALID_INDEX;
    }
  }

  /// \brief Reserves storage for at least \a fromSetSize relation entries.
  void reserve(FromPositionType fromSetSize)
  {
    m_relationsVec.reserve(fromSetSize * relationCardinality());
  }

  void updateSizes()
  {
    m_currentFromSize = m_fromSet->size();
    m_relationsVec.resize(m_currentFromSize * relationCardinality(), INVALID_INDEX);
  }

  /// @}

public:
  /** \brief Direct access to the relation data  */
  RelationVec& data() { return m_relationsVec; }

  /** \brief Direct const access to the relation data  */
  const RelationVec& data() const { return m_relationsVec; }

private:
  inline constexpr FromPositionType relationCardinality() const
  {
    return CardinalityPolicy::size(FromPositionType());
  }

  /**
   * \brief Helper function to expand the relation data storage
   * \param s The requested size
   */
  void expandSizeIfNeeded(FromPositionType s)
  {
    if(s > m_currentFromSize)
    {
      m_currentFromSize = m_fromSet->size();
      m_relationsVec.resize(m_currentFromSize * relationCardinality(), INVALID_INDEX);
    }

    SLIC_ASSERT_MSG(s <= m_currentFromSize,
                    fmt::format("Expanded size {} is larger than relation's 'from' set of {}",
                                s,
                                m_fromSet->size()));
  }

  /**
   * \brief Debug check that an index in the FromSet is not out-of-range
   * \param fromSetIndex An (alleged) index in the FromSet
   */
  inline void verifyPosition(FromPositionType AXOM_DEBUG_PARAM(fromSetIndex)) const
  {
    SLIC_ASSERT_MSG(fromSetIndex >= 0 && fromSetIndex < m_currentFromSize,
                    fmt::format("Index {} out of range [0,{})", fromSetIndex, m_currentFromSize));
  }

private:
  FromSetType* m_fromSet;
  ToSetType* m_toSet;

  RelationVec m_relationsVec;
  IndexType m_currentFromSize {0};
};

/* Checks whether the relation is valid.  */
template <typename PosType, typename ElemType, typename CardinalityPolicy>
bool DynamicConstantRelation<PosType, ElemType, CardinalityPolicy>::isValid(bool verboseOutput) const
{
  fmt::memory_buffer out;

  bool setsAreValid = true;
  bool relationdataIsValid = true;

  // Check if the sets are valid
  const bool isFromSetNull = (m_fromSet == nullptr);
  const bool isToSetNull = (m_toSet == nullptr);

  if(isFromSetNull || isToSetNull)
  {
    if(verboseOutput)
    {
      fmt::format_to(std::back_inserter(out),
                     "\n\t Static relations require both the fromSet and toSet "
                     "to be non-null:"
                     "\t-- fromSet was {} null"
                     "\t-- toSet was {} null",
                     isFromSetNull ? "" : " not ",
                     isToSetNull ? "" : " not ");
    }

    setsAreValid = false;
  }

  // Check the sizes of fromSet matches relationVec
  if(setsAreValid)
  {
    if(m_fromSet->size() != m_currentFromSize)
    {
      if(verboseOutput)
      {
        fmt::format_to(std::back_inserter(out),
                       "\n\t Internal size does not match fromSet size:"
                       "\t-- fromSet size is {}"
                       "\t-- internal size is {}",
                       m_fromSet->size(),
                       m_currentFromSize);
      }
      setsAreValid = false;
    }

    if(m_fromSet->size() * relationCardinality() != (int)m_relationsVec.size())
    {
      if(verboseOutput)
      {
        fmt::format_to(std::back_inserter(out),
                       "\n\t Size of relationVec does not match toSet size:"
                       "\t-- fromSet size is {}"
                       "\t-- m_relationsVec size is {}",
                       m_fromSet->size(),
                       m_relationsVec.size());
      }
      setsAreValid = false;
    }
  }

  // Check if the relation data is valid
  if(setsAreValid)
  {
    // Check that invalid set entry points to invalid relation.
    // Note: the reverse can be valid. ie. valid set entry may have invalid relation entry.
    for(auto pos : m_fromSet->positions())
    {
      if(m_fromSet->at(pos) == FromSetType::INVALID_ENTRY && isValidEntry(pos))
      {
        if(verboseOutput)
        {
          fmt::format_to(std::back_inserter(out),
                         "\n\t* invalid entries in fromSet; has a valid relation at index "
                         "{}, but element not in from set. Values: {}",
                         pos,
                         (*this)[pos]);
        }
        relationdataIsValid = false;
      }
    }

    // Check that all relation indices are in range for m_toSet
    for(auto from_idx : m_fromSet->positions())
    {
      if(m_fromSet->isValidEntry(from_idx))
      {
        for(auto idx = 0; idx < relationCardinality(); ++idx)
        {
          const auto pos = from_idx * relationCardinality() + idx;
          const auto val = m_relationsVec[pos];
          if(val != INVALID_INDEX && !m_toSet->isValidEntry(val))
          {
            if(verboseOutput)
            {
              fmt::format_to(std::back_inserter(out),
                             "\n\t* Relation index out of range or invalid:"
                             "\n\t-- position {} ({}-{}) with value {} needs "
                             "to be in range [0,{}) and index a valid entry",
                             pos,
                             from_idx,
                             idx,
                             val,
                             m_toSet->size());
            }
            relationdataIsValid = false;
          }
        }
      }
    }
  }

  // We are done.  Output the messages if applicable and return
  bool bValid = setsAreValid && relationdataIsValid;

  if(verboseOutput && !bValid)
  {
    SLIC_INFO(fmt::to_string(out));
  }

  return bValid;
}

}  // end namespace axom::slam
