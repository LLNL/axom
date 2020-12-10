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

bool YAMLReader::getValue(const conduit::Node& node, int& value)
{
  // Match LuaReader functionality - narrow from floating-point
  if(node.dtype().is_number())
  {
    value = node.to_int();
    return true;
  }
  return false;
}

bool YAMLReader::getValue(const conduit::Node& node, std::string& value)
{
  if(node.dtype().is_string())
  {
    value = node.as_string();
    return true;
  }
  return false;
}

bool YAMLReader::getValue(const conduit::Node& node, double& value)
{
  // Match LuaReader functionality - promote from integer
  if(node.dtype().is_number())
  {
    value = node.to_float64();
    return true;
  }
  return false;
}
bool YAMLReader::getValue(const conduit::Node& node, bool& value)
{
  if(node.dtype().is_int8())
  {
    value = static_cast<bool>(node.to_int8());
    return true;
  }
  // Boolean literals don't appear to be parsed
  else if(node.dtype().is_string())
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
  return getValue(m_root[id], value);
}

bool YAMLReader::getDouble(const std::string& id, double& value)
{
  return getValue(m_root[id], value);
}

bool YAMLReader::getInt(const std::string& id, int& value)
{
  return getValue(m_root[id], value);
}

bool YAMLReader::getString(const std::string& id, std::string& value)
{
  return getValue(m_root[id], value);
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
  const auto node = m_root[id];
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
  const auto node = m_root[id];
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
  const auto node = m_root[id];
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
  const auto node = m_root[id];
  if(!node.dtype().is_list())
  {
    return false;
  }

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

}  // end namespace inlet
}  // end namespace axom
