// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/**
 * \file Set.hpp
 *
 * \brief Basic API for a Set of entities in a simulation
 *
 */

#include <cstddef>
#include <vector>
#include <type_traits>  // for std::common_type

#include "axom/core/utilities/Utilities.hpp"
#include "axom/slam/Utilities.hpp"

namespace axom::slam
{
/**
 * \class Set
 *
 * \brief Abstract base class for a Set of entities in a simulation
 *
 * Sets can represent mesh vertices, cells, materials or refinement levels.
 * A subset might select boundary vertices or ghost cells. Each position gives
 * access to an element, which need not equal that position. Range-based sets
 * compute their elements, while indirection-backed sets read stored values.
 *
 * This virtual interface provides const positional access and validation.
 * Generic code can instead use SetLike from Concepts.hpp, which does not
 * require inheritance, subsetting or validation methods.
 */
template <typename PosType = slam::DefaultPositionType, typename ElemType = slam::DefaultElementType>
class Set
{
public:
  using PositionType = PosType;
  using ElementType = ElemType;

public:
  // Set () {}
  virtual ~Set() = default;

  /**
   * \brief Random access to the entities of the set
   * \return The value of the element at the given position
   * \pre The position identifies a valid entry in [0, size()).
   */
  [[nodiscard]] virtual ElementType at(PositionType) const = 0;

  /**
   * \brief Get the number of entities in the set
   * \return The number of entities in the set.
   */
  [[nodiscard]] AXOM_HOST_DEVICE virtual PositionType size() const = 0;

  /**
   * \brief Determines if the Set is a Subset of another set.
   * \return true if the set is a subset of another set, otherwise false.
   */
  [[nodiscard]] virtual bool isSubset() const = 0;

  /**
   * \brief Checks whether the set is valid.
   * \return true if the underlying indices are valid, false otherwise.
   */
  [[nodiscard]] virtual bool isValid(bool verboseOutput = false) const = 0;

  /**
   * \brief Return whether size() is zero.
   */
  [[nodiscard]] AXOM_HOST_DEVICE virtual bool empty() const = 0;

#if 0
  /**
   * \brief Returns true if the set contains the given element.
   *
   * Alternatively, we can return the position in the set containing the element,
   * with some value for not containing the element
   */
  virtual bool          contains(const SetElement & elt) const = 0;
#endif

private:
  /// \brief Utility function to verify that the given SetPosition is in a valid range.
  virtual void verifyPosition(PositionType) const = 0;
};

/**
 * \brief General equality operator for two sets.
 * \details Two sets are considered equal if they have the same number of
 * elements and their ordered indices agree.
 */
template <typename P1, typename E1, typename P2, typename E2>
inline bool operator==(const Set<P1, E1>& set1, const Set<P2, E2>& set2)
{
  using PosType = typename std::common_type<P1, P2>::type;
  using ElemType = typename std::common_type<E1, E2>::type;

  PosType const numElts = set1.size();

  // Sets are different if they have a different size
  if(set2.size() != numElts)
  {
    return false;
  }

  // Otherwise, compare the indices element wise
  for(PosType pos = PosType(); pos < numElts; ++pos)
  {
    auto&& e1 = static_cast<ElemType&&>(set1.at(static_cast<P1>(pos)));
    auto&& e2 = static_cast<ElemType&&>(set2.at(static_cast<P2>(pos)));
    if(e1 != e2)
    {
      return false;
    }
  }
  return true;
}
/// \brief Set inequality operator
template <typename P1, typename E1, typename P2, typename E2>
inline bool operator!=(const Set<P1, E1>& set1, const Set<P2, E2>& set2)
{
  return !(set1 == set2);
}

}  // end namespace axom::slam
