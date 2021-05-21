// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_NVTX_INTERFACE_HPP_
#define AXOM_NVTX_INTERFACE_HPP_

#include "axom/core/utilities/nvtx/Macros.hpp"

namespace axom
{
namespace nvtx
{
/*!
 * \brief Predefined set of NVTX colors to use with NVTXRange.
 */
enum class Color : uint32_t
{
  BLACK = 0x00000000,
  GREEN = 0x0000FF00,
  LIME = 0x00BFFF00,
  RED = 0x00FF0000,
  BLUE = 0x000000FF,
  YELLOW = 0x00FFFF00,
  CYAN = 0x0000FFFF,
  MAGENTA = 0x00FF00FF,
  WHITE = 0x00FFFFFF,
  ORANGE = 0x00FFA500,
  PINK = 0x00FF69B4
};

/// \name Glocal Definitions
/// @{

/*!
 * \brief Wildcard used for category
 */
constexpr uint32_t ANY_CATEGORY = 0;

/// @}

/// \name Default Settings
/// @{

/*!
 * \brief Default NVTX color to use. Set to GREEN.
 */
constexpr Color DEFAULT_COLOR = Color::GREEN;

/*!
 * \brief Default NVTX category to use. Set to ANY.
 */
constexpr uint32_t DEFAULT_CATEGORY = ANY_CATEGORY;

/// @}

/// \name NVTX API Functions
/// @{

/*!
 * \brief Sets the color to use with NVTX.
 * \param [in] c the color to use. 
 * 
 * \note If not set, axom::nvtx::DEFAULT_COLOR is used.
 * 
 * \post axom::nvtx::get_color() == c
 */
void set_color(Color c);

/*!
 * \brief Returns the color use by NVTX
 * \return Color the color used by NVTX
 * \see nvtx::Color
 */
Color get_color();

/*!
 * \brief Set the category to use with NVTX.
 * \param [in] category the user-prescribed number for the category to use.
 * 
 * \note If note set, axom::nvtx::DEFAULT_CATEGORY is used.
 * 
 * \post axom::nvtx::get_category() == category
 */
void set_category(uint32_t category);

/*!
 * \brief Get the category used with NVTX.
 * \return category the category used with NVTX.
 */
uint32_t get_category();

/*!
 * \brief Resets the NVTX setting to the defaults.
 * 
 * \post axom::nvtx::get_category() == axom::nvtx::DEFAULT_CATEGORY
 * \post axom::nvtx::get_color() == axom::nvtx::DEFAULT_COLOR
 */
void reset();

/// @}

} /* namespace nvtx */

} /* namespace axom */

#endif /* AXOM_NVTX_INTERFACE_HPP_ */
