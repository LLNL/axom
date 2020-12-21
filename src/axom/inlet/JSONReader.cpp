// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file JSONReader.cpp
 *
 * \brief This file contains the class implementation of the JSONReader.
 *******************************************************************************
 */

#include <fstream>

#include "axom/inlet/JSONReader.hpp"

#include "axom/core/utilities/FileUtilities.hpp"
#include "axom/core/utilities/StringUtilities.hpp"
#include "axom/inlet/inlet_utils.hpp"

#include "fmt/fmt.hpp"
#include "axom/slic.hpp"

// #include "conduit_relay.hpp"

namespace axom
{
namespace inlet
{
JSONReader::JSONReader() { }
bool JSONReader::parseFile(const std::string& filePath)
{
  if(!axom::utilities::filesystem::pathExists(filePath))
  {
    SLIC_WARNING(fmt::format("Inlet: Given JSON input file does not exist: {0}",
                             filePath));
    return false;
  }
  m_root.load(filePath, "json");
  return true;
}

bool JSONReader::parseString(const std::string& JSONString)
{
  if(JSONString.empty())
  {
    SLIC_WARNING("Inlet: Given an empty JSON string to parse.");
    return false;
  }
  m_root.parse(JSONString, "json");
  return true;
}

namespace detail
{
// TODO: allow alternate delimiter at sidre level
const static char SCOPE_DELIMITER = '/';

/*!
 *******************************************************************************
 * \brief Traverses a path starting from a Node
 *
 * \param [in] root The node from which to begin traversal
 * \param [in] id   The path to traverse
 * \note Needed for paths containing integers as a conversion is required
 *******************************************************************************
 */
inline conduit::Node traverseNode(const conduit::Node& root, const std::string& id)
{
  if(root.has_path(id))
  {
    return root[id];
  }

  const conduit::Node* node = &root;
  std::vector<std::string> tokens;
  axom::utilities::string::split(tokens, id, SCOPE_DELIMITER);
  for(const auto& token : tokens)
  {
    // Prefer the string name, but if it doesn't exist, try converting to int
    if(node->has_child(token))
    {
      node = &((*node)[token]);
    }
    else
    {
      auto as_int = checkedConvertToInt(token);
      if(as_int.second && as_int.first < node->number_of_children())
      {
        node = &((*node)[as_int.first]);
      }
      else
      {
        // Bail out and return an empty node to indicate failure
        return {};
      }
    }
  }
  return *node;
}

/*!
 *******************************************************************************
 * \brief Copies a Conduit array into an unordered map, using zero-based array
 * indices as the keys
 *
 * \param [in] array The array to copy from
 * \param [out] map  The map to copy into
 * \note Implementing to allow for widening/narrowing conversions
 *******************************************************************************
 */
template <typename ConduitType, typename MapValueType>
void arrayToMap(const conduit::DataArray<ConduitType>& array,
                std::unordered_map<int, MapValueType>& map)
{
  map.clear();
  for(conduit::index_t i = 0; i < array.number_of_elements(); i++)
  {
    // No begin/end iterators are provided by DataArray
    map[i] = array[i];
  }
}

}  // namespace detail

bool JSONReader::getValue(const conduit::Node& node, int& value)
{
  // Match LuaReader functionality - narrow from floating-point but exclude bool
  if(node.dtype().is_number() && !node.dtype().is_uint8())
  {
    value = node.to_int();
    return true;
  }
  return false;
}

bool JSONReader::getValue(const conduit::Node& node, std::string& value)
{
  if(node.dtype().is_string())
  {
    value = node.as_string();
    return true;
  }
  return false;
}

bool JSONReader::getValue(const conduit::Node& node, double& value)
{
  // Match LuaReader functionality - promote from integer but not bool
  if(node.dtype().is_number() && !node.dtype().is_uint8())
  {
    value = node.to_float64();
    return true;
  }
  return false;
}
bool JSONReader::getValue(const conduit::Node& node, bool& value)
{
  if(node.dtype().is_uint8())
  {
    value = node.as_uint8();
    return true;
  }
  return false;
}

bool JSONReader::getBool(const std::string& id, bool& value)
{
  return getValue(detail::traverseNode(m_root, id), value);
}

bool JSONReader::getDouble(const std::string& id, double& value)
{
  return getValue(detail::traverseNode(m_root, id), value);
}

bool JSONReader::getInt(const std::string& id, int& value)
{
  return getValue(detail::traverseNode(m_root, id), value);
}

bool JSONReader::getString(const std::string& id, std::string& value)
{
  return getValue(detail::traverseNode(m_root, id), value);
}

bool JSONReader::getIntMap(const std::string& id,
                           std::unordered_map<int, int>& values)
{
  return getArray(id, values);
}

bool JSONReader::getDoubleMap(const std::string& id,
                              std::unordered_map<int, double>& values)
{
  return getArray(id, values);
}

bool JSONReader::getBoolMap(const std::string& id,
                            std::unordered_map<int, bool>& values)
{
  return getArray(id, values);
}

bool JSONReader::getStringMap(const std::string& id,
                              std::unordered_map<int, std::string>& values)
{
  return getArray(id, values);
}

bool JSONReader::getIntMap(const std::string& id,
                           std::unordered_map<std::string, int>& values)
{
  return getDictionary(id, values);
}

bool JSONReader::getDoubleMap(const std::string& id,
                              std::unordered_map<std::string, double>& values)
{
  return getDictionary(id, values);
}

bool JSONReader::getBoolMap(const std::string& id,
                            std::unordered_map<std::string, bool>& values)
{
  return getDictionary(id, values);
}

bool JSONReader::getStringMap(const std::string& id,
                              std::unordered_map<std::string, std::string>& values)
{
  return getDictionary(id, values);
}

bool JSONReader::getIndices(const std::string& id, std::vector<int>& indices)
{
  indices.clear();
  const auto node = detail::traverseNode(m_root, id);
  if(!node.dtype().is_list())
  {
    return false;
  }
  indices.resize(node.number_of_children());
  // Arrays in JSON are contiguous so we don't need to query the input file
  std::iota(indices.begin(), indices.end(), 0);
  return true;
}

bool JSONReader::getIndices(const std::string& id,
                            std::vector<std::string>& indices)
{
  indices.clear();
  const auto node = detail::traverseNode(m_root, id);
  if(!node.dtype().is_object())
  {
    return false;
  }
  // FIXME: Update to range-based for loops when Axom begins using Conduit 0.6.0
  auto itr = node.children();
  while(itr.has_next())
  {
    const auto& child = itr.next();
    indices.push_back(child.name());
  }
  return true;
}

FunctionVariant JSONReader::getFunction(const std::string&,
                                        const FunctionType,
                                        const std::vector<FunctionType>&)
{
  SLIC_ERROR("[Inlet] JSON does not support functions");
  return {};
}

template <typename T>
bool JSONReader::getDictionary(const std::string& id,
                               std::unordered_map<std::string, T>& values)
{
  values.clear();
  const auto node = detail::traverseNode(m_root, id);
  if(!node.dtype().is_object())
  {
    return false;
  }

  // FIXME: Update to range-based for loops when Axom begins using Conduit 0.6.0
  auto itr = node.children();
  while(itr.has_next())
  {
    const auto& child = itr.next();
    const auto name = child.name();
    if(!getValue(child, values[name]))
    {
      // The current interface allows for overlapping types, but we need to
      // remove the default-initialized element here if it failed
      values.erase(name);
    }
  }
  return true;
}

template <typename T>
bool JSONReader::getArray(const std::string& id,
                          std::unordered_map<int, T>& values)
{
  values.clear();
  const auto node = detail::traverseNode(m_root, id);
  // Truly primitive (i.e., not string) types are contiguous so we grab the array pointer
  if(node.dtype().number_of_elements() > 1)
  {
    // The template parameter is not enough to know the type of the conduit array
    // as widening/narrowing conversions are supported
    if(node.dtype().is_floating_point())
    {
      detail::arrayToMap(node.as_double_array(), values);
    }
    else if(node.dtype().is_integer())
    {
      detail::arrayToMap(node.as_long_array(), values);
    }
    else
    {
      return false;
    }
  }
  else if(!node.dtype().is_list())
  {
    // Single-element arrays will be just the element itself
    // If it's a single element, we know the index is zero
    if(!getValue(node, values[0]))
    {
      values.erase(0);
      return false;
    }
  }
  // String arrays are not supported natively so the node is directly iterated over
  else
  {
    conduit::index_t index = 0;
    // FIXME: Update to range-based for loops when Axom begins using Conduit 0.6.0
    auto itr = node.children();
    while(itr.has_next())
    {
      const auto& child = itr.next();
      if(!getValue(child, values[index]))
      {
        // The current interface allows for overlapping types, but we need to
        // remove the default-initialized element here if it failed
        values.erase(index);
      }
      index++;
    }
    // Check if nothing was inserted
    if(values.empty())
    {
      return false;
    }
  }
  return true;
}

}  // end namespace inlet
}  // end namespace axom
