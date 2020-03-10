// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file LuaMap.cpp
 *
 * \brief This file contains the class implementation of the LuaMap.
 *******************************************************************************
 */

#include "axom/slim/LuaMap.hpp"

#include "axom/core/utilities/FileUtilities.hpp"
#include "axom/core/utilities/StringUtilities.hpp"

#include "fmt/fmt.hpp"
#include "axom/slic.hpp"

namespace axom
{
namespace slim
{

LuaMap::~LuaMap()
{
  if (m_luaState){
    lua_close(m_luaState);
  }
}


bool LuaMap::parseFile(const std::string& filePath)
{
  if (!axom::utilities::filesystem::pathExists(filePath)) {
    SLIC_WARNING(fmt::format("SLIM: Given Lua input deck does not exist: {0}", filePath));
    return false;
  }

  m_luaState = luaL_newstate();
  if (luaL_loadfile(m_luaState, filePath.c_str()) ||
      lua_pcall(m_luaState, 0, 0, 0)) {
    SLIC_WARNING(fmt::format("SLIM: Given Lua input deck could not be loaded: {0}", filePath));
    m_luaState = nullptr;
    return false;
  }

  return true;
}


bool LuaMap::parseString(const std::string& luaString)
{
  if (luaString.empty()) {
    SLIC_WARNING("SLIM: Given an empty Lua string to parse.");
    return false;
  }

  m_luaState = luaL_newstate();
  if (luaL_loadstring(m_luaState, luaString.c_str()) ||
      lua_pcall(m_luaState, 0, 0, 0)) {
    SLIC_WARNING(fmt::format("SLIM: Given Lua string could not be loaded: {0}", luaString));
    m_luaState = nullptr;
    return false;
  }

  return true;
}


bool LuaMap::findVariable(const std::string& id)
{
  if (!m_luaState){
    SLIC_WARNING("Lua state is not initialized. Call LuaMap::parseString or LuaMap::parseFile first!");
    return false;
  }

  if (axom::utilities::string::startsWith(id, scopeDelimiter)){
    SLIC_WARNING(fmt::format("Variable cannot start with scope delimiter: {0}", id));
    return false;
  }

  if (axom::utilities::string::endsWith(id, scopeDelimiter)){
    SLIC_WARNING(fmt::format("Variable cannot end with scope delimiter: {0}", id));
    return false;
  }

  bool atGlobalScope = true;
  std::vector<std::string> tokens;
  axom::utilities::string::split(tokens, id, scopeDelimiter);

  for (std::string token : tokens) {
    if(atGlobalScope) {
      lua_getglobal(m_luaState, token.c_str());
    } else {
      lua_getfield(m_luaState, -1, token.c_str());
    }
    if(lua_isnil(m_luaState, -1)) {
        SLIC_WARNING(fmt::format("Child variable ('{0}') is not defined in parent variable ('{1}')", token, id));
        return false;
    }
    atGlobalScope = false;
  }

  return true;
}


bool LuaMap::getBool(const std::string& id, bool& value)
{
  if (!findVariable(id)){
    return false;
  }
  value = (bool)lua_toboolean(m_luaState, -1);
  return true;
}


bool LuaMap::getDouble(const std::string& id, double& value)
{
  if (!findVariable(id)){
    return false;
  }
  value = (double)lua_tonumber(m_luaState, -1);
  return true;
}


bool LuaMap::getInt(const std::string& id, int& value)
{
  if (!findVariable(id)){
    return false;
  }
  value = (int)lua_tonumber(m_luaState, -1);
  return true;
}


bool LuaMap::getString(const std::string& id, std::string& value)
{
  if (!findVariable(id)){
    return false;
  }
  value = std::string(lua_tostring(m_luaState, -1));
  return true;
}

} // end namespace slim
} // end namespace axom
