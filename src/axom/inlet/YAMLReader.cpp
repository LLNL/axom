// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file YAMLReader.cpp
 *
 * \brief This file contains the class implementation of the YAMLReader.
 *******************************************************************************
 */

#include <fstream>

#include "axom/inlet/YAMLReader.hpp"

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
YAMLReader::YAMLReader() { }
bool YAMLReader::parseFile(const std::string& filePath)
{
  if(!axom::utilities::filesystem::pathExists(filePath))
  {
    SLIC_WARNING(fmt::format("Inlet: Given YAML input file does not exist: {0}",
                             filePath));
    return false;
  }
  m_root.load(filePath, "yaml");
  return true;
}

bool YAMLReader::parseString(const std::string& YAMLString)
{
  if(YAMLString.empty())
  {
    SLIC_WARNING("Inlet: Given an empty YAML string to parse.");
    return false;
  }
  m_root.parse(YAMLString, "yaml");
  return true;
}

namespace detail
{
template <typename T>
bool hasCompatibleType(const conduit::DataType& type)
{
  return false;
}

template <>
bool hasCompatibleType<int>(const conduit::DataType& type)
{
  // Match LuaReader functionality - allow narrowing from floating-point
  return type.is_number();
}

template <>
bool hasCompatibleType<double>(const conduit::DataType& type)
{
  // Match LuaReader functionality - allow promoting from integer
  return type.is_number();
}

template <>
bool hasCompatibleType<bool>(const conduit::DataType& type)
{
  return type.is_string();
}

template <>
bool hasCompatibleType<std::string>(const conduit::DataType& type)
{
  return type.is_string();
}

template <typename T>
struct yaml_parsed_type
{
  using type = T;
};

template <>
struct yaml_parsed_type<int>
{
  using type = long;
};

template <typename T>
typename yaml_parsed_type<T>::type* getArrayPointer(const conduit::Node&)
{
  return nullptr;
}

template <>
typename yaml_parsed_type<int>::type* getArrayPointer<int>(const conduit::Node& node)
{
  if(!node.dtype().is_number())
  {
    return nullptr;
  }
  // YAML integer literals are parsed as 64-bit ints
  return static_cast<long*>(node.as_long_array().data_ptr());
}

template <>
typename yaml_parsed_type<double>::type* getArrayPointer<double>(
  const conduit::Node& node)
{
  if(!node.dtype().is_number())
  {
    return nullptr;
  }
  return static_cast<double*>(node.as_double_array().data_ptr());
}

// TODO allow alternate delimiter at sidre level
const static char SCOPE_DELIMITER = '/';
conduit::Node getNodeChild(const conduit::Node& root, const std::string& id)
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
}  // namespace detail

bool YAMLReader::getValue(const conduit::Node& node, int& value)
{
  // Match LuaReader functionality - narrow from floating-point
  if(detail::hasCompatibleType<int>(node.dtype()))
  {
    value = node.to_int();
    return true;
  }
  return false;
}

bool YAMLReader::getValue(const conduit::Node& node, std::string& value)
{
  if(detail::hasCompatibleType<std::string>(node.dtype()))
  {
    value = node.as_string();
    return true;
  }
  return false;
}

bool YAMLReader::getValue(const conduit::Node& node, double& value)
{
  // Match LuaReader functionality - promote from integer
  if(detail::hasCompatibleType<double>(node.dtype()))
  {
    value = node.to_float64();
    return true;
  }
  return false;
}
bool YAMLReader::getValue(const conduit::Node& node, bool& value)
{
  // Boolean literals don't appear to be parsed as such - they are strings
  if(node.dtype().is_string())
  {
    // YAML 1.2 spec, section 10.3.2
    std::string as_str = node.as_string();
    std::transform(as_str.begin(),
                   as_str.end(),
                   as_str.begin(),
                   [](const unsigned char c) { return std::tolower(c); });
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
  return false;
}

bool YAMLReader::getBool(const std::string& id, bool& value)
{
  return getValue(detail::getNodeChild(m_root, id), value);
}

bool YAMLReader::getDouble(const std::string& id, double& value)
{
  return getValue(detail::getNodeChild(m_root, id), value);
}

bool YAMLReader::getInt(const std::string& id, int& value)
{
  return getValue(detail::getNodeChild(m_root, id), value);
}

bool YAMLReader::getString(const std::string& id, std::string& value)
{
  return getValue(detail::getNodeChild(m_root, id), value);
}

bool YAMLReader::getIntMap(const std::string& id,
                           std::unordered_map<int, int>& values)
{
  return getArray(id, values);
}

bool YAMLReader::getDoubleMap(const std::string& id,
                              std::unordered_map<int, double>& values)
{
  return getArray(id, values);
}

bool YAMLReader::getBoolMap(const std::string& id,
                            std::unordered_map<int, bool>& values)
{
  return getArray(id, values);
}

bool YAMLReader::getStringMap(const std::string& id,
                              std::unordered_map<int, std::string>& values)
{
  return getArray(id, values);
}

bool YAMLReader::getIntMap(const std::string& id,
                           std::unordered_map<std::string, int>& values)
{
  return getDictionary(id, values);
}

bool YAMLReader::getDoubleMap(const std::string& id,
                              std::unordered_map<std::string, double>& values)
{
  return getDictionary(id, values);
}

bool YAMLReader::getBoolMap(const std::string& id,
                            std::unordered_map<std::string, bool>& values)
{
  return getDictionary(id, values);
}

bool YAMLReader::getStringMap(const std::string& id,
                              std::unordered_map<std::string, std::string>& values)
{
  return getDictionary(id, values);
}

bool YAMLReader::getIndices(const std::string& id, std::vector<int>& indices)
{
  indices.clear();
  const auto node = detail::getNodeChild(m_root, id);
  if(!node.dtype().is_list())
  {
    return false;
  }
  indices.resize(node.number_of_children());
  std::iota(indices.begin(), indices.end(), 0);
  return true;
}

bool YAMLReader::getIndices(const std::string& id,
                            std::vector<std::string>& indices)
{
  indices.clear();
  const auto node = detail::getNodeChild(m_root, id);
  if(!node.dtype().is_object())
  {
    return false;
  }
  for(const auto& child : node.children())
  {
    indices.push_back(child.name());
  }
  return true;
}

template <typename T>
bool YAMLReader::getDictionary(const std::string& id,
                               std::unordered_map<std::string, T>& values)
{
  values.clear();
  const auto node = detail::getNodeChild(m_root, id);
  if(!node.dtype().is_object())
  {
    return false;
  }

  for(const auto& child : node.children())
  {
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
bool YAMLReader::getArray(const std::string& id,
                          std::unordered_map<int, T>& values)
{
  values.clear();
  const auto node = detail::getNodeChild(m_root, id);
  // Truly primitive (i.e., not string) types are contiguous so we grab the array pointer
  if(node.dtype().number_of_elements() > 1)
  {
    // A layer of indirection is needed here - cannot instantiate a conduit::DataArray
    if(auto data_ptr = detail::getArrayPointer<T>(node))
    {
      for(conduit::index_t i = 0; i < node.dtype().number_of_elements(); i++)
      {
        // No begin/end iterators provided for conduit::DataArray
        values[i] = data_ptr[i];
      }
      return true;
    }
    return false;
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
    return true;
  }
  // String arrays are not supported natively so they cacn be directly iterated over
  else
  {
    conduit::index_t index = 0;
    for(const auto& child : node.children())
    {
      if(!getValue(child, values[index]))
      {
        // The current interface allows for overlapping types, but we need to
        // remove the default-initialized element here if it failed
        values.erase(index);
      }
      index++;
    }
    return true;
  }
}

}  // end namespace inlet
}  // end namespace axom
