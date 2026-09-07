// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/**
 * \file DynamicMap.hpp
 * \brief Map values that can grow with a dynamic set.
 */

#include <vector>
#include <sstream>

#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"
#include "axom/slic/interface/slic.hpp"

#include "axom/slam/DynamicSet.hpp"
#include "axom/slam/Map.hpp"

namespace axom::slam
{
/**
 * \class DynamicMap
 * \brief Store one value per set position, with storage that can grow.
 *
 * An entry is valid when the corresponding set entry is valid. The map owns its
 * std::vector of values and borrows the set. Mutable operator[] can grow the
 * value buffer, while flatValue() requires an existing valid entry. Neither
 * operation adds elements to the set, which must outlive the map.
 */
template <typename SetT, typename DataT>
class DynamicMap
{
public:
  using SetType = SetT;
  /// The complete set bound by set().
  using MappedSetType = SetType;
  using DataType = DataT;

  using PositionType = typename SetType::PositionType;
  using SetElement = typename SetType::ElementType;
  using ValueType = DataType&;
  using ConstValueType = const DataType&;

  using OrderedMap = std::vector<DataType>;

public:
  /** \brief Default constructor   */
  DynamicMap() : m_set(nullptr) { }

  /**
   * \brief Constructor from a set pointer
   *
   * \param theSet The set, which must outlive the map.
   *
   * Allocates theSet->size() value-initialized entries for a non-null set.
   */
  DynamicMap(SetType* theSet) : m_set(theSet)
  {
    if(m_set != nullptr)
    {
      m_data.resize(m_set->size());
    }
  }

  /**
   * \brief Constructor from a set pointer
   *
   * \param theSet The set, which must outlive the map.
   * \param defaultValue Initial value of every entry.
   *
   * The map will be allocated with \a theSet->size() entries.
   * Each entry will have value \a defaultValue
   */
  DynamicMap(SetType* theSet, DataType defaultValue) : m_set(theSet)
  {
    if(m_set != nullptr)
    {
      m_data.resize(m_set->size(), defaultValue);
    }
  }

  ~DynamicMap() { }

public:
  /// \brief Returns a pointer to the map's underlying set
  const SetType* set() const { return m_set; }

  /// \name DynamicMap individual access functions
  /// @{
  ///

  /// \brief Return the value at set index \a setIndex
  const DataType& operator[](PositionType setIndex) const
  {
    verifyPosition(setIndex);
    return m_data[setIndex];
  }

  /// \brief Access the scalar component of an existing, valid entry.
  /// \pre component == 0. Unlike mutable operator[], this does not grow the map.
  DataType& flatValue(PositionType pos, PositionType AXOM_DEBUG_PARAM(component) = 0)
  {
    SLIC_ASSERT(component == 0);
    verifyPosition(pos);
    return m_data[pos];
  }

  /// \overload
  const DataType& flatValue(PositionType pos, PositionType AXOM_DEBUG_PARAM(component) = 0) const
  {
    SLIC_ASSERT(component == 0);
    return (*this)[pos];
  }

  /// \brief Return the set element associated with a valid entry.
  SetElement index(PositionType pos) const
  {
    verifyPosition(pos);
    return set()->at(pos);
  }

  /// \brief DynamicMap has one scalar component per entry.
  PositionType numComp() const { return 1; }

  /// @}

  /// \brief Access to underlying data
  OrderedMap& data() { return m_data; }

  /// \brief Const access to underlying data
  const OrderedMap& data() const { return m_data; }

  /// \brief Reserves storage for at least \a s entries.
  void reserve(PositionType s) { m_data.reserve(s); }

  /// \name DynamicMap cardinality functions
  /// @{

  /** \brief Returns the size of map's set */
  PositionType size() const { return static_cast<PositionType>(m_data.size()); }

  /**
   * \brief Return the number of valid entries
   *
   * An entry at a given index is considered valid if corresponding set element is valid.
   */
  PositionType numberOfValidEntries() const
  {
    return (m_set != nullptr) ? m_set->numberOfValidEntries() : 0;
  }

  /// @}

  /// \name DynamicMap validity check functions
  /// @{

  bool isValidEntry(PositionType pos) const
  {
    return (m_set != nullptr) ? m_set->isValidEntry(pos) : false;
  }

  /// \brief Predicate to check if this DynamicMap instance is valid
  [[nodiscard]] bool isValid(bool verboseOutput = false) const;

  /// @}

private:
  /// \brief Debug check that the index is not out of range
  inline void verifyPosition(PositionType AXOM_DEBUG_PARAM(setIndex)) const
  {
    SLIC_ASSERT_MSG(
      setIndex >= 0 && setIndex < size(),
      "Attempted to access entry " << setIndex << " but map's set has size " << m_data.size());

    SLIC_ASSERT_MSG(isValidEntry(setIndex), "Attempted to access invalid set entry " << setIndex);
  }

public:
  /// \name Functions that modify the map's cardinality
  /// @{

  /**
   * \brief Get the value at position \a position
   *
   * \note Increases the map size if position is out of range
   */
  DataType& operator[](PositionType position)
  {
    if(size() < position + 1)
    {
      resize(position + 1);
    }

    return m_data[position];
  }

  /**
   * \brief Insert \a value into \a position in the map
   *
   * \note Increases the map size if position is out of range
   */
  void insert(PositionType position, DataType value) { operator[](position) = value; }

  /**
   * \brief Resizes the map to have at least \a s positions
   * \param s The minimum necessary capacity for resizing
   * \pre s >= 0
   */
  void resize(PositionType s)
  {
    // Note (KW): For this to be a valid DynamicMap operation,
    // we would need to also increase the size of the map's set!

    // Note (KW): Do we want this to shrink the size of the map
    // when s < size() ?

    SLIC_ASSERT_MSG(s >= 0, "Attempted to resize vector with a negative size " << s);

    m_data.resize(s);
  }

  /// @}

private:
  SetType* m_set;
  OrderedMap m_data;
};

template <typename SetT, typename DataT>
bool DynamicMap<SetT, DataT>::isValid(bool verboseOutput) const
{
  bool bValid = true;

  std::stringstream errStr;

  if(m_set == nullptr)
  {
    if(!m_data.empty())
    {
      if(verboseOutput)
      {
        errStr << "\n\t* the underlying set was never provided,"
               << " but its associated data is not empty"
               << " , data has size " << m_data.size();
      }
      bValid = false;
    }
  }
  else
  {
    // Check the data array and set data have equal size
    if(static_cast<PositionType>(m_data.size()) != m_set->size())
    {
      if(verboseOutput)
      {
        errStr << "\n\t* the underlying set and its associated mapped data"
               << " have different sizes, underlying set has size " << m_set->size()
               << " , data has size " << m_data.size();
        ;
      }

      bValid = false;
    }
  }

  if(verboseOutput)
  {
    if(bValid)
    {
      SLIC_DEBUG("Map was valid.");
    }
    else
    {
      SLIC_DEBUG("Map was not valid. " << errStr.str());
    }
  }

  return bValid;
}

}  // end namespace axom::slam
