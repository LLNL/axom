// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file LuaReader.cpp
 *
 * \brief This file contains the class implementation of the LuaReader.
 *******************************************************************************
 */

#include "axom/inlet/LuaReader.hpp"

#include "axom/core/utilities/FileUtilities.hpp"
#include "axom/core/utilities/StringUtilities.hpp"

#include "fmt/fmt.hpp"
#include "axom/slic.hpp"

namespace axom
{
namespace inlet
{

LuaReader::~LuaReader()
{
  if (m_luaState)
  {
    lua_close(m_luaState);
  }
}


bool LuaReader::parseFile(const std::string& filePath)
{
  if (!axom::utilities::filesystem::pathExists(filePath))
  {
    SLIC_WARNING(fmt::format("Inlet: Given Lua input deck does not exist: {0}",
                             filePath));
    return false;
  }

  if (!m_luaState)
    m_luaState = luaL_newstate();
  
  if (luaL_loadfile(m_luaState, filePath.c_str()) ||
      lua_pcall(m_luaState, 0, 0, 0))
  {
    SLIC_WARNING(fmt::format(
                   "Inlet: Given Lua input deck could not be loaded: {0}",
                   filePath));
    m_luaState = nullptr;
    return false;
  }

  return true;
}


bool LuaReader::parseString(const std::string& luaString)
{
  if (luaString.empty())
  {
    SLIC_WARNING("Inlet: Given an empty Lua string to parse.");
    return false;
  }

  if (!m_luaState)
    m_luaState = luaL_newstate();

  if (luaL_loadstring(m_luaState, luaString.c_str()) ||
      lua_pcall(m_luaState, 0, 0, 0))
  {
    SLIC_WARNING(fmt::format("Inlet: Given Lua string could not be loaded: {0}",
                             luaString));
    m_luaState = nullptr;
    return false;
  }

  return true;
}

// TODO allow alternate delimiter at sidre level
#define SCOPE_DELIMITER '/'

bool LuaReader::findVariable(const std::string& id)
{
  if (!m_luaState)
  {
    SLIC_WARNING(
      "Lua state is not initialized. Call LuaReader::parseString or LuaReader::parseFile first!");
    return false;
  }

  std::string temp_id = id;
  //TODO: support multiple roots?
  if (axom::utilities::string::startsWith(temp_id, SCOPE_DELIMITER))
  {
    temp_id.erase(0, 1);
  }

  if (axom::utilities::string::endsWith(id, SCOPE_DELIMITER))
  {
    SLIC_WARNING(fmt::format("Variable cannot end with scope delimiter: {0}",
                             id));
    return false;
  }

  bool atGlobalScope = true;
  std::vector<std::string> tokens;
  axom::utilities::string::split(tokens, temp_id, SCOPE_DELIMITER);

  // Clear the lua stack because we always call with fully qualified names
  lua_settop(m_luaState, 0);
  for (std::string token : tokens)
  {
    if(atGlobalScope)
    {
      lua_getglobal(m_luaState, token.c_str());
    }
    else
    {
      lua_getfield(m_luaState, -1, token.c_str());
    }
    if(lua_isnil(m_luaState, -1))
    {
      // variable not found
      return false;
    }
    atGlobalScope = false;
  }

  return true;
}


bool LuaReader::getBool(const std::string& id, bool& value)
{
  if (!findVariable(id))
  {
    return false;
  }
  value = (bool)lua_toboolean(m_luaState, -1);
  return true;
}


bool LuaReader::getDouble(const std::string& id, double& value)
{
  if (!findVariable(id))
  {
    return false;
  }
  value = (double)lua_tonumber(m_luaState, -1);
  return true;
}


bool LuaReader::getInt(const std::string& id, int& value)
{
  if (!findVariable(id))
  {
    return false;
  }
  value = (int)lua_tonumber(m_luaState, -1);
  return true;
}


bool LuaReader::getString(const std::string& id, std::string& value)
{
  if (!findVariable(id))
  {
    return false;
  }
  value = std::string(lua_tostring(m_luaState, -1));
  return true;
}

  void LuaReader::releaseLuaState()
  {
    m_luaState = nullptr;
  }
  
} // end namespace inlet
} // end namespace axom
