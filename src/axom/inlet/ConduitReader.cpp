// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file ConduitReader.cpp
 *
 * \brief This file contains the class implementation of the ConduitReader.
 *******************************************************************************
 */

#include "axom/inlet/ConduitReader.hpp"

#include <fstream>

#include "axom/core/utilities/FileUtilities.hpp"
#include "axom/core/utilities/StringUtilities.hpp"
#include "axom/inlet/inlet_utils.hpp"

#include "fmt/fmt.hpp"
#include "axom/slic.hpp"

namespace axom
{
namespace inlet
{
ConduitReader::ConduitReader(const std::string& protocol) : m_protocol(protocol)
{
  SLIC_ERROR_IF((protocol != "yaml") && (protocol != "json"),
                fmt::format("Inlet: Only 'json' and 'yaml' protocols are "
                            "supported by ConduitReader, got: {0}",
                            protocol));
}

bool ConduitReader::parseFile(const std::string& filePath)
{
  if(!axom::utilities::filesystem::pathExists(filePath))
  {
    SLIC_WARNING(
      fmt::format("Inlet: Given input file does not exist: {0}", filePath));
    return false;
  }
  bool success = true;
  // Temporarily enable exceptions for Conduit errors so they can be caught and displayed
  sidre::DataStore::setConduitDefaultMessageHandlers();
  try
  {
    m_root.load(filePath, m_protocol);
  }
  catch(const conduit::Error& e)
  {
    SLIC_WARNING(
      fmt::format("[Inlet]: Failed to parse {0}:\n{1}", m_protocol, e.message()));
    success = false;
  }
  sidre::DataStore::setConduitSLICMessageHandlers();
  return success;
}

bool ConduitReader::parseString(const std::string& stringToRead)
{
  if(stringToRead.empty())
  {
    SLIC_WARNING("Inlet: Given an empty string to parse.");
    return false;
  }
  bool success = true;
  sidre::DataStore::setConduitDefaultMessageHandlers();
  try
  {
    m_root.parse(stringToRead, m_protocol);
  }
  catch(const conduit::Error& e)
  {
    SLIC_WARNING(
      fmt::format("[Inlet]: Failed to parse {0}:\n{1}", m_protocol, e.message()));
    success = false;
  }
  sidre::DataStore::setConduitSLICMessageHandlers();
  return success;
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
 * \param [in] id   The path to traverse as an Inlet path
 * \note Needed for paths containing integers as a conversion is required
 * \note Inlet paths differ from Conduit paths in that they can contain integers,
 * e.g., "/path/to/7/foo" would correspond to node["path"]["to"][7]["foo"].
 * The traversal prefers the string-valued name but will attempt to convert to
 * integer if no such child exists.
 *******************************************************************************
 */
conduit::Node traverseNode(const conduit::Node& root, const std::string& id)
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

bool ConduitReader::getValue(const conduit::Node& node, int& value)
{
  // Match LuaReader functionality - narrow from floating-point but exclude bool
  if(node.dtype().is_number() && !node.dtype().is_uint8())
  {
    value = node.to_int();
    return true;
  }
  return false;
}

bool ConduitReader::getValue(const conduit::Node& node, std::string& value)
{
  if(node.dtype().is_string())
  {
    value = node.as_string();
    return true;
  }
  return false;
}

bool ConduitReader::getValue(const conduit::Node& node, double& value)
{
  // Match LuaReader functionality - promote from integer but not bool
  if(node.dtype().is_number() && !node.dtype().is_uint8())
  {
    value = node.to_double();
    return true;
  }
  return false;
}

bool ConduitReader::getValue(const conduit::Node& node, bool& value)
{
  // Boolean literals don't appear to be parsed as such - they are strings
  if((m_protocol == "yaml") && node.dtype().is_string())
  {
    std::string as_str = node.as_string();
    // YAML 1.2 spec, section 10.3.2
    // FIXME: Converting the string to lowercase is not strictly correct, it
    // allows for things like tRue and falsE
    utilities::string::toLower(as_str);
    if(as_str == "true")
    {
      value = true;
      return true;
    }
    else if(as_str == "false")
    {
      value = false;
      return true;
    }
  }
  else if((m_protocol == "json") && node.dtype().is_uint8())
  {
    value = node.as_uint8();
    return true;
  }
  return false;
}

bool ConduitReader::getBool(const std::string& id, bool& value)
{
  return getValue(detail::traverseNode(m_root, id), value);
}

bool ConduitReader::getDouble(const std::string& id, double& value)
{
  return getValue(detail::traverseNode(m_root, id), value);
}

bool ConduitReader::getInt(const std::string& id, int& value)
{
  return getValue(detail::traverseNode(m_root, id), value);
}

bool ConduitReader::getString(const std::string& id, std::string& value)
{
  return getValue(detail::traverseNode(m_root, id), value);
}

bool ConduitReader::getIntMap(const std::string& id,
                              std::unordered_map<int, int>& values)
{
  return getArray(id, values);
}

bool ConduitReader::getDoubleMap(const std::string& id,
                                 std::unordered_map<int, double>& values)
{
  return getArray(id, values);
}

bool ConduitReader::getBoolMap(const std::string& id,
                               std::unordered_map<int, bool>& values)
{
  return getArray(id, values);
}

bool ConduitReader::getStringMap(const std::string& id,
                                 std::unordered_map<int, std::string>& values)
{
  return getArray(id, values);
}

bool ConduitReader::getIntMap(const std::string& id,
                              std::unordered_map<VariantKey, int>& values)
{
  return getDictionary(id, values);
}

bool ConduitReader::getDoubleMap(const std::string& id,
                                 std::unordered_map<VariantKey, double>& values)
{
  return getDictionary(id, values);
}

bool ConduitReader::getBoolMap(const std::string& id,
                               std::unordered_map<VariantKey, bool>& values)
{
  return getDictionary(id, values);
}

bool ConduitReader::getStringMap(const std::string& id,
                                 std::unordered_map<VariantKey, std::string>& values)
{
  return getDictionary(id, values);
}

bool ConduitReader::getIndices(const std::string& id, std::vector<int>& indices)
{
  indices.clear();
  const auto node = detail::traverseNode(m_root, id);
  int num_elements = node.number_of_children();
  // Primitive arrays do not count as lists
  if(!node.dtype().is_list())
  {
    num_elements = node.dtype().number_of_elements();
  }

  indices.resize(num_elements);
  // Arrays in YAML/JSON are contiguous so we don't need to query the input file
  std::iota(indices.begin(), indices.end(), 0);
  return true;
}

bool ConduitReader::getIndices(const std::string& id,
                               std::vector<VariantKey>& indices)
{
  indices.clear();
  const auto node = detail::traverseNode(m_root, id);
  if(!node.dtype().is_object())
  {
    // If it's not an object, try integer indexing
    std::vector<int> int_indices;
    if(getIndices(id, int_indices))
    {
      for(const int idx : int_indices)
      {
        indices.emplace_back(idx);
      }
      return true;
    }
    return false;
  }
  for(const auto& child : node.children())
  {
    indices.push_back(child.name());
  }
  return true;
}

FunctionVariant ConduitReader::getFunction(const std::string&,
                                           const FunctionType,
                                           const std::vector<FunctionType>&)
{
  SLIC_ERROR("[Inlet] Conduit YAML/JSON does not support functions");
  return {};
}

template <typename T>
bool ConduitReader::getDictionary(const std::string& id,
                                  std::unordered_map<VariantKey, T>& values)
{
  values.clear();
  const auto node = detail::traverseNode(m_root, id);
  if(!node.dtype().is_object())
  {
    return false;
  }

  for(const auto& child : node.children())
  {
    const auto name = child.name();

    T value;
    // Inlet allows for heterogenous containers, so a failure here is "normal"
    if(getValue(child, value))
    {
      values[name] = value;
    }
  }
  return true;
}

template <typename T>
bool ConduitReader::getArray(const std::string& id,
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
    T value;
    if(getValue(node, value))
    {
      values[0] = value;
    }
    else
    {
      return false;
    }
  }
  // String arrays are not supported natively so the node is directly iterated over
  else
  {
    conduit::index_t index = 0;
    for(const auto& child : node.children())
    {
      T value;
      // Inlet allows for heterogenous containers, so a failure here is "normal"
      if(getValue(child, value))
      {
        values[index] = value;
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
