// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
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
  LuaReader() { m_lua.open_libraries(sol::lib::base); }

  /*!
   *****************************************************************************
   * \brief Parses the given input file.
   *
   * This performs any setup work and parses the given input file.
   * It is required that this is called before using the Reader and overrides
   * any Lua state that was previously there.
   *
   * \param [in] filePath The Input file to be read
   *
   * \return true if the input file was able to be parsed
   *****************************************************************************
   */
  bool parseFile(const std::string& filePath);

  /*!
   *****************************************************************************
   * \brief Parses the given Lua string.
   *
   * This performs any setup work and parses the given Lua string.
   * It is required that this is called before using the Reader and overrides
   * any Lua state that was previously there.
   *
   * \param [in] luaString The Input file to be read
   *
   * \return true if the string was able to be parsed
   *****************************************************************************
   */
  bool parseString(const std::string& luaString);

  /*!
   *****************************************************************************
   * \brief Return a boolean out of the input file
   *
   * This performs any necessary retrieval and mapping from the given identifier
   * to what is in the input file.
   *
   * \param [in]  id    The identifier to the bool that will be retrieved
   * \param [out] value The value of the bool that was retrieved
   *
   * \return true if the variable was able to be retrieved from the file
   *****************************************************************************
   */
  bool getBool(const std::string& id, bool& value);

  /*!
   *****************************************************************************
   * \brief Return a double out of the input file
   *
   * This performs any necessary retrieval and mapping from the given identifier
   * to what is in the input file.
   *
   * \param [in]  id    The identifier to the double that will be retrieved
   * \param [out] value The value of the double that was retrieved
   *
   * \return true if the variable was able to be retrieved from the file
   *****************************************************************************
   */
  bool getDouble(const std::string& id, double& value);

  /*!
   *****************************************************************************
   * \brief Return a int out of the input file
   *
   * This performs any necessary retrieval and mapping from the given identifier
   * to what is in the input file.
   *
   * \param [in]  id    The identifier to the int that will be retrieved from
   *the file
   * \param [out] value The value of the int that was retrieved from the file
   *
   * \return true if the variable was able to be retrieved from the file
   *****************************************************************************
   */
  bool getInt(const std::string& id, int& value);

  /*!
   *****************************************************************************
   * \brief Return a string out of the input file
   *
   * This performs any necessary retrieval and mapping from the given identifier
   * to what is in the input file.
   *
   * \param [in]  id    The identifier to the string that will be retrieved
   * \param [out] value The value of the string that was retrieved
   *
   * \return true if the variable was able to be retrieved from the file
   *****************************************************************************
   */
  bool getString(const std::string& id, std::string& value);

  /*!
   *****************************************************************************
   * \brief Get an index-integer mapping for the given Lua array
   *
   * This performs any necessary retrieval and mapping from the given identifier
   * to what is in the input file.
   *
   * \param [in]  id    The identifier to the string that will be retrieved
   * \param [out] map The values of the ints that were retrieved
   *
   * \return true if the array was able to be retrieved from the file
   *****************************************************************************
   */
  bool getIntMap(const std::string& id, std::unordered_map<int, int>& values);
  bool getIntMap(const std::string& id,
                 std::unordered_map<std::string, int>& values);

  /*!
   *****************************************************************************
   * \brief Get an index-double mapping for the given Lua array
   *
   * This performs any necessary retrieval and mapping from the given identifier
   * to what is in the input file.
   *
   * \param [in]  id    The identifier to the string that will be retrieved
   * \param [out] map The values of the doubles that were retrieved
   *
   * \return true if the array was able to be retrieved from the file
   *****************************************************************************
   */
  bool getDoubleMap(const std::string& id,
                    std::unordered_map<int, double>& values);
  bool getDoubleMap(const std::string& id,
                    std::unordered_map<std::string, double>& values);

  /*!
   *****************************************************************************
   * \brief Get an index-bool mapping for the given Lua array
   *
   * This performs any necessary retrieval and mapping from the given identifier
   * to what is in the input file.
   *
   * \param [in]  id    The identifier to the string that will be retrieved
   * \param [out] map The values of the bools that were retrieved
   *
   * \return true if the array was able to be retrieved from the file
   *****************************************************************************
   */
  bool getBoolMap(const std::string& id, std::unordered_map<int, bool>& values);
  bool getBoolMap(const std::string& id,
                  std::unordered_map<std::string, bool>& values);

  /*!
   *****************************************************************************
   * \brief Get an index-string mapping for the given Lua array
   *
   * This performs any necessary retrieval and mapping from the given identifier
   * to what is in the input file.
   *
   * \param [in]  id    The identifier to the string that will be retrieved
   * \param [out] map The values of the strings that were retrieved
   *
   * \return true if the array was able to be retrieved from the file
   *****************************************************************************
   */
  bool getStringMap(const std::string& id,
                    std::unordered_map<int, std::string>& values);
  bool getStringMap(const std::string& id,
                    std::unordered_map<std::string, std::string>& values);

  /*!
   *****************************************************************************
   * \brief Get the list of indices for an container
   *
   * \param [in]  id    The identifier to the container that will be retrieved
   * \param [out] map The values of the indices that were retrieved
   *
   * \return true if the indices were able to be retrieved from the file
   *****************************************************************************
   */
  bool getIndices(const std::string& id, std::vector<int>& indices);
  bool getIndices(const std::string& id, std::vector<std::string>& indices);

  /*!
   *****************************************************************************
   * \brief Get a function from the input deck
   *
   * \param [in]  id    The identifier to the function that will be retrieved
   * \param [in]  ret_type    The return type of the function
   * \param [in]  arg_types    The argument types of the function
   *
   * \return The function, compares false if not found
   *****************************************************************************
   */
  FunctionVariant getFunction(const std::string& id,
                              const FunctionTag ret_type,
                              const std::vector<FunctionTag>& arg_types);

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
  bool getValue(const std::string& id, T& value);

  // Expect this to be called for only Inlet-supported types.
  template <typename Key, typename Val>
  bool getMap(const std::string& id,
              std::unordered_map<Key, Val>& values,
              sol::type type);

  template <typename T>
  bool getIndicesInternal(const std::string& id, std::vector<T>& indices);

  /*!
   *****************************************************************************
   * \brief Obtains the Lua table reached by successive indexing through the
   * container of keys described by a pair of iterators
   * 
   * \note For a set of keys {key1, key2, key3, ...}, this function
   * is equivalent to
   * \code{.cpp}
   * table = m_lua[key1][key2][key3][...];
   * \endcode
   * 
   * \param [in] begin Iterator to the beginning of the container of keys
   * \param [in] end Iterator to one-past-then-end of the container
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
};

}  // end namespace inlet
}  // end namespace axom

#endif
