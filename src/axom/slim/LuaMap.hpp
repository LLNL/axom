// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file LuaMap.hpp
 *
 * \brief This file contains the class definition of the LuaMap.
 *******************************************************************************
 */

#ifndef LUAMAP_HPP
#define LUAMAP_HPP

#include "axom/slim/Map.hpp"

extern "C" {
  #include "lua.h"
  #include "lualib.h"
  #include "lauxlib.h"
}

namespace axom
{
namespace slim
{

/*!
 *******************************************************************************
 * \class LuaMap
 *
 * \brief A Map that is able to read and map the variables from a Lua input deck.
 *
 * \see Map
 *******************************************************************************
 */
class LuaMap : public Map
{
public:
  /*!
   *****************************************************************************
   * \brief Destructor for the LuaMap.
   *
   * This performs any cleanup work the Map needs to do before going
   * away.
   *****************************************************************************
   */
  ~LuaMap();

  /*!
   *****************************************************************************
   * \brief Parses the given input deck.
   *
   * This performs any setup work and parses the given input deck.
   * It is required that this is called before using the Map and overrides
   * any Lua state that was previously there.
   *
   * \param [in] filePath The Input deck to be read
   *
   * \return true if the input deck was able to be parsed
   *****************************************************************************
   */
  bool parseFile(const std::string& filePath);

  /*!
   *****************************************************************************
   * \brief Parses the given Lua string.
   *
   * This performs any setup work and parses the given Lua string.
   * It is required that this is called before using the Map and overrides
   * any Lua state that was previously there.
   *
   * \param [in] luaString The Input deck to be read
   *
   * \return true if the string was able to be parsed
   *****************************************************************************
   */
  bool parseString(const std::string& luaString);

  /*!
   *****************************************************************************
   * \brief Return a boolean out of the input deck
   *
   * This performs any necessary retrieval and mapping from the given identifier
   * to what is in the input deck.
   *
   * \param [in]  id    The identifier to the bool that will be retrieved from the deck
   * \param [out] value The value of the bool that was retrieved from the deck
   *
   * \return true if the variable was able to be retrieved from the deck
   *****************************************************************************
   */
  bool getBool(const std::string& id, bool& value);

  /*!
   *****************************************************************************
   * \brief Return a double out of the input deck
   *
   * This performs any necessary retrieval and mapping from the given identifier
   * to what is in the input deck.
   *
   * \param [in]  id    The identifier to the double that will be retrieved from the deck
   * \param [out] value The value of the double that was retrieved from the deck
   *
   * \return true if the variable was able to be retrieved from the deck
   *****************************************************************************
   */
  bool getDouble(const std::string& id, double& value);

  /*!
   *****************************************************************************
   * \brief Return a int out of the input deck
   *
   * This performs any necessary retrieval and mapping from the given identifier
   * to what is in the input deck.
   *
   * \param [in]  id    The identifier to the int that will be retrieved from the deck
   * \param [out] value The value of the int that was retrieved from the deck
   *
   * \return true if the variable was able to be retrieved from the deck
   *****************************************************************************
   */
  bool getInt(const std::string& id, int& value);

  /*!
   *****************************************************************************
   * \brief Return a string out of the input deck
   *
   * This performs any necessary retrieval and mapping from the given identifier
   * to what is in the input deck.
   *
   * \param [in]  id    The identifier to the string that will be retrieved from the deck
   * \param [out] value The value of the string that was retrieved from the deck
   *
   * \return true if the variable was able to be retrieved from the deck
   *****************************************************************************
   */
  bool getString(const std::string& id, std::string& value);

private:

  /*!
   *****************************************************************************
   * \brief Move the Lua state to the given Lua variable id
   *
   * This performs any necessary retrieval and mapping from the given identifier
   * to what is in the input deck.
   *
   * \param [in] id The identifier to the bool that will be retrieved from the deck
   *
   * \return true if Lua variable id was found and Lua state was moved to variable
   *****************************************************************************
   */
  bool findVariable(const std::string& id);

  lua_State* m_luaState;
};

} // end namespace slim
} // end namespace axom

#endif