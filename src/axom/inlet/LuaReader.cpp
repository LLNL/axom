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

#include <fstream>

#include "axom/inlet/LuaReader.hpp"

#include "axom/core/utilities/FileUtilities.hpp"
#include "axom/core/utilities/StringUtilities.hpp"
#include "axom/inlet/inlet_utils.hpp"

#include "fmt/fmt.hpp"
#include "axom/slic.hpp"

namespace axom
{
namespace inlet
{
bool LuaReader::parseFile(const std::string& filePath)
{
  if(!axom::utilities::filesystem::pathExists(filePath))
  {
    SLIC_WARNING(
      fmt::format("Inlet: Given Lua input file does not exist: {0}", filePath));
    return false;
  }

  auto script = m_lua.script_file(filePath);
  if(!script.valid())
  {
    SLIC_WARNING(
      fmt::format("Inlet: Given Lua input deck is invalid: {0}", filePath));
  }
  return script.valid();
}

bool LuaReader::parseString(const std::string& luaString)
{
  if(luaString.empty())
  {
    SLIC_WARNING("Inlet: Given an empty Lua string to parse.");
    return false;
  }
  m_lua.script(luaString);
  return true;
}

// TODO allow alternate delimiter at sidre level
#define SCOPE_DELIMITER '/'

bool LuaReader::getBool(const std::string& id, bool& value)
{
  return getValue(id, value);
}

bool LuaReader::getDouble(const std::string& id, double& value)
{
  return getValue(id, value);
}

bool LuaReader::getInt(const std::string& id, int& value)
{
  return getValue(id, value);
}

bool LuaReader::getString(const std::string& id, std::string& value)
{
  return getValue(id, value);
}

bool LuaReader::getIntMap(const std::string& id,
                          std::unordered_map<int, int>& values)
{
  return getMap(id, values, sol::type::number);
}

bool LuaReader::getDoubleMap(const std::string& id,
                             std::unordered_map<int, double>& values)
{
  return getMap(id, values, sol::type::number);
}

bool LuaReader::getBoolMap(const std::string& id,
                           std::unordered_map<int, bool>& values)
{
  return getMap(id, values, sol::type::boolean);
}

bool LuaReader::getStringMap(const std::string& id,
                             std::unordered_map<int, std::string>& values)
{
  return getMap(id, values, sol::type::string);
}

bool LuaReader::getArrayIndices(const std::string& id, std::vector<int>& indices)
{
  indices.clear();
  std::vector<std::string> tokens;
  axom::utilities::string::split(tokens, id, SCOPE_DELIMITER);

  if(tokens.empty() || !m_lua[tokens[0]].valid())
  {
    return false;
  }
  sol::table t = m_lua[tokens[0]];
  for(size_t i = 1; i < tokens.size(); i++)
  {
    // Use the C versions to avoid the exceptions
    // thrown by std::stoi on conversion failure
    char* ptr;
    auto as_int = strtol(tokens[i].c_str(), &ptr, 10);
    if((!*ptr) && t[as_int].valid())
    {
      t = t[as_int];
    }
    else if(t[tokens[i]].valid())
    {
      t = t[tokens[i]];
    }
    else
    {
      return false;
    }
  }

  for(const auto& entry : t)
  {
    indices.push_back(entry.first.as<int>());
  }
  return true;
}

template <typename T>
bool LuaReader::getValue(const std::string& id, T& value)
{
  std::vector<std::string> tokens;
  axom::utilities::string::split(tokens, id, SCOPE_DELIMITER);
  if(tokens.size() == 1 && m_lua[id].valid())
  {
    value = m_lua[id];
    return true;
  }

  if(!m_lua[tokens[0]].valid())
  {
    return false;
  }
  sol::table t = m_lua[tokens[0]];
  for(size_t i = 1; i < tokens.size() - 1; i++)
  {
    // Use the C versions to avoid the exceptions
    // thrown by std::stoi on conversion failure
    char* ptr;
    auto as_int = strtol(tokens[i].c_str(), &ptr, 10);
    if((!*ptr) && t[as_int].valid())
    {
      t = t[as_int];
    }
    else if(t[tokens[i]].valid())
    {
      t = t[tokens[i]];
    }
    else
    {
      return false;
    }
  }
  if(t[tokens.back()].valid())
  {
    value = t[tokens.back()];
    return true;
  }

  return false;
}

template <typename T>
bool LuaReader::getMap(const std::string& id,
                       std::unordered_map<int, T>& values,
                       sol::type type)
{
  values.clear();
  std::vector<std::string> tokens;
  axom::utilities::string::split(tokens, id, SCOPE_DELIMITER);

  if(tokens.empty() || !m_lua[tokens[0]].valid())
  {
    return false;
  }
  sol::table t = m_lua[tokens[0]];
  for(size_t i = 1; i < tokens.size(); i++)
  {
    // Use the C versions to avoid the exceptions
    // thrown by std::stoi on conversion failure
    char* ptr;
    auto as_int = strtol(tokens[i].c_str(), &ptr, 10);
    if((!*ptr) && t[as_int].valid())
    {
      t = t[as_int];
    }
    else if(t[tokens[i]].valid())
    {
      t = t[tokens[i]];
    }
    else
    {
      return false;
    }
  }

  for(const auto& entry : t)
  {
    // Gets only indexed items in the table.
    if(entry.first.get_type() == sol::type::number &&
       entry.second.get_type() == type)
    {
      values[entry.first.as<int>()] = entry.second.as<T>();
    }
  }
  return true;
}

}  // end namespace inlet
}  // end namespace axom
