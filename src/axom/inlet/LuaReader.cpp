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
      fmt::format("Inlet: Given Lua input file is invalid: {0}", filePath));
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

template <typename Iter>
bool LuaReader::traverseToTable(Iter begin, Iter end, sol::table& table)
{
  // Nothing to traverse
  if(begin == end)
  {
    return true;
  }

  if(!m_lua[*begin].valid())
  {
    return false;
  }

  table = m_lua[*begin];  // Use the first one to index into the global lua state
  ++begin;

  // Then use the remaining keys to walk down to the requested table
  for(auto curr = begin; curr != end; ++curr)
  {
    auto key = *curr;
    // Use the C versions to avoid the exceptions
    // thrown by std::stoi on conversion failure
    // FIXME: Switch to std::from_chars when C++17 is available
    char* ptr;
    auto as_int = strtol(key.c_str(), &ptr, 10);
    if((!*ptr) && table[as_int].valid())
    {
      table = table[as_int];
    }
    else if(table[key].valid())
    {
      table = table[key];
    }
    else
    {
      return false;
    }
  }
  return true;
}

bool LuaReader::getArrayIndices(const std::string& id, std::vector<int>& indices)
{
  std::vector<std::string> tokens;
  axom::utilities::string::split(tokens, id, SCOPE_DELIMITER);

  sol::table t;

  if(tokens.empty() || !traverseToTable(tokens.begin(), tokens.end(), t))
  {
    return false;
  }

  indices.clear();

  // std::transform ends up being messier here
  for(const auto& entry : t)
  {
    indices.push_back(entry.first.as<int>());
  }
  return true;
}

sol::protected_function LuaReader::getFunctionInternal(const std::string& id)
{
  std::vector<std::string> tokens;
  axom::utilities::string::split(tokens, id, SCOPE_DELIMITER);

  sol::protected_function lua_func;

  if(tokens.size() == 1)
  {
    if(m_lua[tokens[0]].valid())
    {
      lua_func = m_lua[tokens[0]];
    }
  }
  else
  {
    sol::table t;
    // Don't traverse through the last token as it doesn't contain a table
    if(traverseToTable(tokens.begin(), tokens.end() - 1, t) &&
       t[tokens.back()].valid())
    {
      lua_func = t[tokens.back()];
    }
  }
  return lua_func;
}

namespace detail
{
template <typename Arg>
sol::protected_function_result callWith(sol::protected_function func,
                                        const Arg& arg)
{
  auto tentative_result = func(arg);
  SLIC_ERROR_IF(!tentative_result.valid(), "[Inlet] Lua function call failed");
  return tentative_result;
}

sol::protected_function_result callWith(sol::protected_function func,
                                        const primal::Vector2D& vec)
{
  auto tentative_result = func(vec[0], vec[1]);
  SLIC_ERROR_IF(!tentative_result.valid(), "[Inlet] Lua function call failed");
  return tentative_result;
}

sol::protected_function_result callWith(sol::protected_function func,
                                        const primal::Vector3D& vec)
{
  auto tentative_result = func(vec[0], vec[1], vec[2]);
  SLIC_ERROR_IF(!tentative_result.valid(), "[Inlet] Lua function call failed");
  return tentative_result;
}

template <typename Ret>
Ret extractResult(sol::protected_function_result&& res)
{
  return res;
}

template <>
primal::Vector2D extractResult<primal::Vector2D>(sol::protected_function_result&& res)
{
  primal::Vector2D result;
  std::tuple<int, int> tup = res;
  result[0] = std::get<0>(tup);
  result[1] = std::get<1>(tup);
  return result;
}

template <>
primal::Vector3D extractResult<primal::Vector3D>(sol::protected_function_result&& res)
{
  primal::Vector3D result;
  std::tuple<int, int, int> tup = res;
  result[0] = std::get<0>(tup);
  result[1] = std::get<1>(tup);
  result[2] = std::get<2>(tup);
  return result;
}

template <typename Ret, typename Arg>
std::function<Ret(typename inlet_function_arg_type<Arg>::type)> buildStdFunction(
  sol::protected_function&& func)
{
  return [func](typename inlet_function_arg_type<Arg>::type arg) {
    return extractResult<Ret>(callWith(std::move(func), arg));
  };
}

template <InletFunctionType Ret>
InletFunctionWrapper bindArgType(sol::protected_function&& func,
                                 const InletFunctionType arg_type)
{
  switch(arg_type)
  {
  case InletFunctionType::Vec2D:
    return buildStdFunction<typename inlet_function_type<Ret>::type, primal::Vector2D>(
      std::move(func));
  case InletFunctionType::Vec3D:
    return buildStdFunction<typename inlet_function_type<Ret>::type, primal::Vector3D>(
      std::move(func));
  case InletFunctionType::Double:
    return buildStdFunction<typename inlet_function_type<Ret>::type, double>(
      std::move(func));
  default:
    SLIC_ERROR("[Inlet] Unexpected function argument type");
  }
  return {};  // Never reached but needed as errors do not imply control flow as with exceptions
}

}  // end namespace detail

InletFunctionWrapper LuaReader::getFunction(const std::string& id,
                                            const InletFunctionType ret_type,
                                            const InletFunctionType arg_type)
{
  auto lua_func = getFunctionInternal(id);
  SLIC_ERROR_IF(
    !lua_func,
    "[Inlet] Function did not exist at the specified path in the input file");
  switch(ret_type)
  {
  case InletFunctionType::Vec2D:
    return detail::bindArgType<InletFunctionType::Vec2D>(std::move(lua_func),
                                                         arg_type);
  case InletFunctionType::Vec3D:
    return detail::bindArgType<InletFunctionType::Vec3D>(std::move(lua_func),
                                                         arg_type);
  case InletFunctionType::Double:
    return detail::bindArgType<InletFunctionType::Double>(std::move(lua_func),
                                                          arg_type);
  default:
    SLIC_ERROR("[Inlet] Unexpected function return type");
  }
  return {};  // Never reached but needed as errors do not imply control flow as with exceptions
}

template <typename T>
bool LuaReader::getValue(const std::string& id, T& value)
{
  std::vector<std::string> tokens;
  axom::utilities::string::split(tokens, id, SCOPE_DELIMITER);

  if(tokens.size() == 1)
  {
    if(m_lua[tokens[0]].valid())
    {
      value = m_lua[tokens[0]];
      return true;
    }
    return false;
  }

  sol::table t;
  // Don't traverse through the last token as it doesn't contain a table
  if(!traverseToTable(tokens.begin(), tokens.end() - 1, t))
  {
    return false;
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

  sol::table t;
  if(tokens.empty() || !traverseToTable(tokens.begin(), tokens.end(), t))
  {
    return false;
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
