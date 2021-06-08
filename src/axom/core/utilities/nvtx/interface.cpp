// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/utilities/nvtx/interface.hpp"

namespace axom
{
namespace nvtx
{
/*!
 * \brief Internal data-structure to hold the NVTX settings.
 */
static struct settings_t
{
  uint32_t color;    /*!< the associated color with an NVTX Range */
  uint32_t category; /*!< the associated category with an NVTX Range */

  settings_t()
    : color(static_cast<uint32_t>(DEFAULT_COLOR))
    , category(DEFAULT_CATEGORY)
  { }

} Settings;

//------------------------------------------------------------------------------
// NVTX interface implementation
//------------------------------------------------------------------------------

void set_color(Color c) { Settings.color = static_cast<uint32_t>(c); }

//------------------------------------------------------------------------------
Color get_color() { return static_cast<Color>(Settings.color); }

//------------------------------------------------------------------------------
void set_category(uint32_t category) { Settings.category = category; }

//------------------------------------------------------------------------------
uint32_t get_category() { return Settings.category; }

//------------------------------------------------------------------------------
void reset()
{
  set_color(DEFAULT_COLOR);
  set_category(DEFAULT_CATEGORY);
}

} /* namespace nvtx */

} /* namespace axom */
