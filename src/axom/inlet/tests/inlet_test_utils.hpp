// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef INLET_TEST_UTILS_HPP
#define INLET_TEST_UTILS_HPP

#include <vector>
#include <string>

#include "gtest/gtest.h"

#ifdef AXOM_USE_SOL
  #include "axom/inlet/LuaReader.hpp"
#endif

namespace axom::inlet::detail
{
/*!
 *******************************************************************************
 * \brief Converts a Lua string to a string accepted by the reader of specified
 * type
 * \param [in] luaString The Lua code representing the contents of an input file
 * \tparam InletReader The inlet::Reader implementation intended to receive
 * the result of this function
 *******************************************************************************
 */
template <typename InletReader>
inline std::string fromLuaTo(const std::string& luaString)
{
  return luaString;
}

#ifdef AXOM_USE_SOL
using ReaderTypes = ::testing::Types<axom::inlet::LuaReader>;
#else
using ReaderTypes = ::testing::Types<>;
#endif

}  // namespace axom::inlet::detail

#endif
