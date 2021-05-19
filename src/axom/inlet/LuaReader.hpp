// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file LuaReader.hpp
 *
 * \brief This file contains the class definition of the LuaReader.
 *******************************************************************************
 */

#ifndef INLET_LUAMAP_HPP
#define INLET_LUAMAP_HPP

#include "axom/inlet/Reader.hpp"
#include <sol/sol.hpp>

extern "C" {
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"
}

namespace axom
{
namespace inlet
{
/*!
 *******************************************************************************
 * \class LuaReader
 *
 * \brief A Reader that is able to read variables from a Lua file.
 *
 * \see Reader
 *******************************************************************************
 */
class LuaReader : public Reader
{
public:
  LuaReader();

  bool parseFile(const std::string& filePath) override;

  bool parseString(const std::string& luaString) override;

  ReaderResult getBool(const std::string& id, bool& value) override;

  ReaderResult getDouble(const std::string& id, double& value) override;

  ReaderResult getInt(const std::string& id, int& value) override;

  ReaderResult getString(const std::string& id, std::string& value) override;

  ReaderResult getIntMap(const std::string& id,
                         std::unordered_map<int, int>& values) override;
  ReaderResult getIntMap(const std::string& id,
                         std::unordered_map<VariantKey, int>& values) override;

  ReaderResult getDoubleMap(const std::string& id,
                            std::unordered_map<int, double>& values) override;
  ReaderResult getDoubleMap(const std::string& id,
                            std::unordered_map<VariantKey, double>& values) override;

  ReaderResult getBoolMap(const std::string& id,
                          std::unordered_map<int, bool>& values) override;
  ReaderResult getBoolMap(const std::string& id,
                          std::unordered_map<VariantKey, bool>& values) override;

  ReaderResult getStringMap(const std::string& id,
                            std::unordered_map<int, std::string>& values) override;
  ReaderResult getStringMap(
    const std::string& id,
    std::unordered_map<VariantKey, std::string>& values) override;

  ReaderResult getIndices(const std::string& id,
                          std::vector<int>& indices) override;
  ReaderResult getIndices(const std::string& id,
                          std::vector<VariantKey>& indices) override;

  FunctionVariant getFunction(const std::string& id,
                              const FunctionTag ret_type,
                              const std::vector<FunctionTag>& arg_types) override;

  std::vector<std::string> getAllNames() override;

  /*!
   *****************************************************************************
   * \brief The base index for arrays in Lua
   *****************************************************************************
   */
  static const int baseIndex = 1;

  /*!
   *****************************************************************************
   * \brief Returns the Sol Lua state
   *
   * This allows the user to access functionality that was not provided by Inlet.
   *
   * \return Reference to the Sol Lua state
   *****************************************************************************
   */
  sol::state& solState() { return m_lua; }

private:
  // Expect this to be called for only Inlet-supported types.
  template <typename T>
  ReaderResult getValue(const std::string& id, T& value);

  // Expect this to be called for only Inlet-supported types.
  template <typename Key, typename Val>
  ReaderResult getMap(const std::string& id,
                      std::unordered_map<Key, Val>& values,
                      sol::type type);

  template <typename T>
  ReaderResult getIndicesInternal(const std::string& id, std::vector<T>& indices);

  /*!
   *****************************************************************************
   * \brief Obtains the Lua table reached by successive indexing through the
   * range of keys described by a pair of iterators
   * 
   * \note For a set of keys {key1, key2, key3, ...}, this function
   * is equivalent to
   * \code{.cpp}
   * table = m_lua[key1][key2][key3][...];
   * \endcode
   * 
   * \param [in] begin Iterator to the beginning of the range of keys
   * \param [in] end Iterator to one-past-the-end of the range
   * \param [out] t The table to traverse
   * 
   * \return Whether the traversal was successful
   *****************************************************************************
   */
  template <typename Iter>
  bool traverseToTable(Iter begin, Iter end, sol::table& table);

  /*!
   *****************************************************************************
   * \brief Traverses the Lua state to retrieve a sol function object
   *
   * \param [in]  id    The identifier to the function that will be retrieved
   *
   * \return The function, compares false if not found
   *****************************************************************************
   */
  sol::protected_function getFunctionInternal(const std::string& id);

  sol::state m_lua;
  // The elements in the global table preloaded by Sol/Lua, these are ignored
  // to ensure that name retrieval only includes user-provided paths
  std::vector<std::string> m_preloaded_globals;
};

}  // end namespace inlet
}  // end namespace axom

#endif
