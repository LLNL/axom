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
LuaReader::LuaReader()
{
  m_lua.open_libraries(sol::lib::base);
  auto vec_type = m_lua.new_usertype<primal::Vector3D>(
    "Vec3D",  // Name of the class in Lua
    // Add make_vector as a constructor to enable "new Vec3D(x,y,z)"
    // Use lambdas for 2D and "default" cases - default arguments cannot be
    // propagated automatically
    "new",
    sol::factories(
      primal::Vector3D::make_vector,
      [](double x, double y) {
        return primal::Vector3D {x, y};
      },
      [] { return primal::Vector3D {}; }),
    // Add vector addition operation
    sol::meta_function::addition,
    sol::resolve<primal::Vector3D(const primal::Vector3D&, const primal::Vector3D&)>(
      primal::operator+),
    // Needs to be resolved in the same way as operator+
    sol::meta_function::unary_minus,
    sol::resolve<primal::Vector3D(const primal::Vector3D&)>(primal::operator-),
    // To allow both "directions" of a scalar multiplication, the overloads
    // have to be manually specified + resolved
    sol::meta_function::multiplication,
    sol::overload(
      sol::resolve<primal::Vector3D(const primal::Vector3D&, const double)>(
        primal::operator*),
      sol::resolve<primal::Vector3D(const double, const primal::Vector3D&)>(
        primal::operator*)),
    // Separate functions from get/set via index - subtract 1 as lua is 1-indexed
    sol::meta_function::index,
    [](const primal::Vector3D& vec, const int key) { return vec[key - 1]; },
    // A lambda is used here as the set-via-returned reference is insufficient
    sol::meta_function::new_index,
    [](primal::Vector3D& vec, const int key, const double value) {
      vec[key - 1] = value;
    },
    // Set up the mathematical operations by name
    "norm",
    &primal::Vector3D::norm,
    "squared_norm",
    &primal::Vector3D::squared_norm,
    "unitVector",
    &primal::Vector3D::unitVector,
    "dot",
    &primal::Vector3D::dot,
    // Implemented as a member function like dot products are, for simplicity
    "cross",
    // Needs to be resolved as it is an overloaded static method
    sol::resolve<primal::Vector3D(const primal::Vector3D&, const primal::Vector3D&)>(
      primal::Vector3D::cross_product));
}

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

// A set of pure functions for handling the conversion of Lua functions to C++
// callables
namespace detail
{
// Passes through everything except a vector, which is expanded into
// three separate arguments
template <typename Arg>
Arg&& lua_identity(Arg&& arg)
{
  return std::forward<Arg>(arg);
}

std::tuple<double, double, double> lua_identity(const primal::Vector3D& vec)
{
  return std::make_tuple(vec[0], vec[1], vec[2]);
}

/*!
 *****************************************************************************
 * \brief Templated function for calling a sol function
 *
 * \param [in] func The sol function of unknown concrete type
 * \tparam Args The argument types of the function
 *
 * \return A checkable version of the function's result
 *****************************************************************************
 */
template <typename... Args>
sol::protected_function_result callWith(const sol::protected_function& func,
                                        Args&&... args)
{
  auto tentative_result = func(lua_identity(std::forward<Args>(args))...);
  SLIC_ERROR_IF(
    !tentative_result.valid(),
    "[Inlet] Lua function call failed, argument types possibly incorrect");
  return tentative_result;
}

/*!
 *****************************************************************************
 * \brief Templated function for extracting a concrete type from a sol function
 * result, used to allow for returning nonprimitive types, specifically, vectors
 *
 * \param [in] res The sol result of unknown concrete type
 * \tparam Ret The return type of the function
 *
 * \return The function's result
 *****************************************************************************
 */
template <typename Ret>
Ret extractResult(sol::protected_function_result&& res)
{
  sol::optional<Ret> option = res;
  SLIC_ERROR_IF(
    !option,
    "[Inlet] Lua function call failed, return types possibly incorrect");
  return option.value();
}

template <>
primal::Vector3D extractResult<primal::Vector3D>(sol::protected_function_result&& res)
{
  auto tup = extractResult<std::tuple<double, double, double>>(std::move(res));
  return {std::get<0>(tup), std::get<1>(tup), std::get<2>(tup)};
}

/*!
 *****************************************************************************
 * \brief Creates a std::function given a Lua function and template parameters
 * corresponding to the function signature
 *
 * \param [in] func The sol object containing the lua function of unknown signature
 * \tparam Ret The return type of the function
 * \tparam Args... The argument types of the function
 *
 * \return A std::function that wraps the lua function
 * 
 * \note This is needed as a layer of indirection for bindArgType so it can
 * properly deduce the constructor call
 *****************************************************************************
 */
template <typename Ret, typename... Args>
std::function<Ret(typename detail::inlet_function_arg_type<Args>::type...)>
buildStdFunction(sol::protected_function&& func)
{
  // Generalized lambda capture needed to move into lambda
  return [func(std::move(func))](
           typename detail::inlet_function_arg_type<Args>::type... args) {
    return extractResult<Ret>(callWith(func, args...));
  };
}

/*!
 *****************************************************************************
 * \brief Adds argument types to a parameter pack based on the contents
 * of a std::vector of type tags
 *
 * \param [in] func The sol object containing the lua function of unknown signature
 * \param [in] arg_types The vector of argument types
 * 
 * \tparam I The number of arguments processed, or "stack size", used to mitigate
 * infinite compile-time recursion
 * \tparam Ret The function's return type
 * \tparam Args... The function's current arguments (already processed), remaining
 * arguments are in the arg_types vector
 *
 * \return A callable wrapper
 *****************************************************************************
 */
template <std::size_t I, typename Ret, typename... Args>
typename std::enable_if<(I > MAX_NUM_ARGS), FunctionVariant>::type bindArgType(
  sol::protected_function&&,
  const std::vector<FunctionType>&)
{
  SLIC_ERROR("[Inlet] Maximum number of function arguments exceeded: " << I);
  return {};
}

template <std::size_t I, typename Ret, typename... Args>
typename std::enable_if<I <= MAX_NUM_ARGS, FunctionVariant>::type bindArgType(
  sol::protected_function&& func,
  const std::vector<FunctionType>& arg_types)
{
  if(arg_types.size() == I)
  {
    return buildStdFunction<Ret, Args...>(std::move(func));
  }
  else
  {
    switch(arg_types[I])
    {
    case FunctionType::Vec3D:
      return bindArgType<I + 1, Ret, Args..., primal::Vector3D>(std::move(func),
                                                                arg_types);
    case FunctionType::Double:
      return bindArgType<I + 1, Ret, Args..., double>(std::move(func), arg_types);
    default:
      SLIC_ERROR("[Inlet] Unexpected function argument type");
    }
  }
  return {};  // Never reached but needed as errors do not imply control flow as with exceptions
}

/*!
 *****************************************************************************
 * \brief Performs a type-checked access to a Lua table
 *
 * \param [in]  proxy The sol::proxy object to retrieve from
 * \param [out] val The value to write to, if it is of the correct type
 *
 * \return true if the value was retrieved from the lua state
 *****************************************************************************
 */
template <typename Proxy, typename Value>
bool checkedGet(const Proxy& proxy, Value& val)
{
  sol::optional<Value> option = proxy;
  if(option)
  {
    val = option.value();
    return true;
  }
  return false;
}

}  // end namespace detail

FunctionVariant LuaReader::getFunction(const std::string& id,
                                       const FunctionType ret_type,
                                       const std::vector<FunctionType>& arg_types)
{
  auto lua_func = getFunctionInternal(id);
  if(lua_func)
  {
    switch(ret_type)
    {
    case FunctionType::Vec3D:
      return detail::bindArgType<0u, primal::Vector3D>(std::move(lua_func),
                                                       arg_types);
    case FunctionType::Double:
      return detail::bindArgType<0u, double>(std::move(lua_func), arg_types);
    default:
      SLIC_ERROR("[Inlet] Unexpected function return type");
    }
  }
  return {};  // Return an empty function to indicate that the function was not found
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
      return detail::checkedGet(m_lua[tokens[0]], value);
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
    return detail::checkedGet(t[tokens.back()], value);
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
      detail::checkedGet(m_lua[tokens[0]], lua_func);
    }
  }
  else
  {
    sol::table t;
    // Don't traverse through the last token as it doesn't contain a table
    if(traverseToTable(tokens.begin(), tokens.end() - 1, t) &&
       t[tokens.back()].valid())
    {
      detail::checkedGet(t[tokens.back()], lua_func);
    }
  }
  return lua_func;
}

}  // end namespace inlet
}  // end namespace axom
