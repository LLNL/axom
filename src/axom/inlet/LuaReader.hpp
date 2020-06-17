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

#ifndef INLET_LUAMAP_HPP
#define INLET_LUAMAP_HPP

#include "axom/inlet/Reader.hpp"

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
 * \brief A Reader that is able to read variables from a Lua deck.
 *
 * \see Reader
 *******************************************************************************
 */
class LuaReader : public Reader
{
public:
  /*!
   *****************************************************************************
   * \brief Destructor for the LuaReader.
   *
   * This performs any cleanup work the Reader needs to do before going
   * away.
   *****************************************************************************
   */
  ~LuaReader();

  /*!
   *****************************************************************************
   * \brief Parses the given input deck.
   *
   * This performs any setup work and parses the given input deck.
   * It is required that this is called before using the Reader and overrides
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
   * It is required that this is called before using the Reader and overrides
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
   * \param [in]  id    The identifier to the bool that will be retrieved
   * \param [out] value The value of the bool that was retrieved
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
   * \param [in]  id    The identifier to the double that will be retrieved
   * \param [out] value The value of the double that was retrieved
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
   * \param [in]  id    The identifier to the int that will be retrieved from
   *the deck
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
   * \param [in]  id    The identifier to the string that will be retrieved
   * \param [out] value The value of the string that was retrieved
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
   * This performs any necessary retrieval and mapping from the given fully
   * qualified Sidre path name to what is in the input deck.
   *
   * \param [in] id The path name to the bool that will be retrieved
   *
   * \return true if Lua variable was found and Lua state was moved to
   * variable
   *****************************************************************************
   */
  bool findVariable(const std::string& id);

  lua_State* m_luaState;
};

} // end namespace inlet
} // end namespace axom

#endif
