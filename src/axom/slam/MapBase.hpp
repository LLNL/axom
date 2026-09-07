// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/**
 * \file MapBase.hpp
 *
 * \brief The virtual size and validation interface for maps.
 *
 */

#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"

#include "axom/slam/Set.hpp"

namespace axom::slam
{
/**
 * \class   MapBase
 *
 * \brief Query map size and validity through a common base pointer.
 *
 * Value access belongs to the concrete map type. MapLike in Concepts.hpp
 * describes the operations used by generic algorithms without requiring this base.
 * \see     Map
 *
 */

template <typename SetPositionType = slam::DefaultPositionType>
class MapBase
{
public:
  using PositionType = SetPositionType;

public:
  AXOM_HOST_DEVICE
  virtual ~MapBase() { };

  /**
   * \brief Get the number of entities in the set used by this map
   * \return The number of entities in the set used in the map.
   */
  [[nodiscard]] AXOM_HOST_DEVICE virtual PositionType size() const = 0;

  /**
   * \brief Checks whether the Map is valid.
   * \return   True if valid, false otherwise.
   */
  [[nodiscard]] virtual bool isValid(bool verboseOutput) const = 0;

private:
  /**
   * \brief Verifies that the provided PositionType is in a valid range.
   */
  virtual void verifyPosition(PositionType) const = 0;
};

}  // end namespace axom::slam
