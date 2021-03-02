// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file JSONSchemaWriter.cpp
 *
 * \brief This file contains the class implementation of the JSONSchemaWriter.
 *******************************************************************************
 */

#include "axom/inlet/JSONSchemaWriter.hpp"

#include <iostream>

#include "axom/slic.hpp"
#include "axom/inlet/Container.hpp"

namespace axom
{
namespace inlet
{
namespace detail
{
sidre::TypeID recordRange(const sidre::View* view, conduit::Node& node)
{
  auto type = view->getTypeID();
  if(type == sidre::INT_ID)
  {
    const int* range = view->getData();
    node["minimum"] = range[0];
    node["maximum"] = range[1];
  }
  else
  {
    const double* range = view->getData();
    node["minimum"] = range[0];
    node["maximum"] = range[1];
  }

  return type;
}

sidre::TypeID recordEnum(const sidre::View* view, conduit::Node& node)
{
  sidre::TypeID type = view->getTypeID();
  auto size = view->getNumElements();
  if(type == sidre::INT_ID)
  {
    const int* values = view->getData();
    for(int i = 0; i < size; i++)
    {
      node["enum"].append() = values[i];
    }
  }
  else
  {
    const double* values = view->getData();
    for(int i = 0; i < size; i++)
    {
      node["enum"].append() = values[i];
    }
  }
  return type;
}

sidre::TypeID recordEnum(const sidre::Group* group, conduit::Node& node)
{
  auto idx = group->getFirstValidViewIndex();
  while(axom::sidre::indexIsValid(idx))
  {
    node["enum"].append() = std::string(group->getView(idx)->getString());
    idx = group->getNextValidViewIndex(idx);
  }
  return sidre::CHAR8_STR_ID;
}

void recordFieldSchema(const Field& field, conduit::Node& node)
{
  const static auto jsonTypeNames = []() {
    std::unordered_map<sidre::TypeID, std::string> result;
    result[sidre::TypeID::INT_ID] = "integer";
    result[sidre::TypeID::CHAR8_STR_ID] = "string";
    result[sidre::TypeID::FLOAT64_ID] = "number";
    result[sidre::TypeID::INT8_ID] = "boolean";
    return result;
  }();

  // The type information isn't guaranteed to be anywhere, so we try the
  // value, the valid values, the range, etc
  axom::sidre::TypeID type = sidre::NO_TYPE_ID;

  const auto sidreGroup = field.sidreGroup();

  if(sidreGroup->hasView("description"))
  {
    node["description"] =
      std::string(sidreGroup->getView("description")->getString());
  }

  if(sidreGroup->hasView("range"))
  {
    type = recordRange(sidreGroup->getView("range"), node);
  }
  else if(sidreGroup->hasView("validValues"))
  {
    type = recordEnum(sidreGroup->getView("validValues"), node);
  }
  else if(sidreGroup->hasGroup("validStringValues"))
  {
    type = recordEnum(sidreGroup->getGroup("validStringValues"), node);
  }

  if(sidreGroup->hasView("defaultValue"))
  {
    auto default_view = sidreGroup->getView("defaultValue");
    type = default_view->getTypeID();
    switch(type)
    {
    case sidre::INT_ID:
      node["default"] = static_cast<int>(default_view->getData());
      break;
    case sidre::FLOAT64_ID:
      node["default"] = static_cast<double>(default_view->getData());
      break;
    case sidre::CHAR8_STR_ID:
      node["default"] = std::string(default_view->getString());
      break;
    case sidre::INT8_ID:
      node["default"] = static_cast<int8>(default_view->getData());
      break;
    default:
      break;
    }
  }

  if(sidreGroup->hasView("value"))
  {
    type = sidreGroup->getView("value")->getTypeID();
  }

  // If we couldn't figure out what type it is, we can't enforce it in the schema
  if(type != sidre::NO_TYPE_ID)
  {
    node["type"] = jsonTypeNames.at(type);
  }
}

}  // namespace detail

JSONSchemaWriter::JSONSchemaWriter(const std::string& fileName)
  : m_fileName(fileName)
{
  m_schemaRoot["$schema"] = "https://json-schema.org/draft/2020-12/schema#";
  m_schemaRoot["title"] = "Input File Options";
  m_schemaRoot["description"] = "The schema for the input file";
  m_schemaRoot["type"] = "object";
}

void JSONSchemaWriter::documentContainer(const Container& container)
{
  const auto sidreGroup = container.sidreGroup();
  // Ignore everything before the root - intended for situations when the Inlet
  // root is not the root of the datastore
  if(!m_rootPathInitialized)
  {
    m_rootPath = sidreGroup->getPath();
    m_rootPathInitialized = true;
  }

  // FIXME: Use path class for all of this logic - all this is doing is inserting
  // a 'properties' in after each token
  std::vector<std::string> tokens;
  const std::string filteredPathName =
    removePrefix(m_rootPath, sidreGroup->getPathName());
  utilities::string::split(tokens, filteredPathName, '/');
  std::string containerPath =
    fmt::format("properties/{}", fmt::join(tokens, "/properties/"));
  auto& containerNode =
    container.name().empty() ? m_schemaRoot : m_schemaRoot[containerPath];
  // if (!tokens.empty())
  // {
  //   tokens.pop_back(); // Remove the basename so we can add it later
  // }
  // auto& containerNode = m_schemaRoot;
  // for (const auto& token : tokens)
  // {
  //   containerNode = containerNode[token]["properties"];
  // }
  // if (!container.name().empty())
  // {
  //   // Re-add the basename to get the "root" of the container
  //   containerNode = containerNode[sidreGroup->getName()];
  // }
  // const std::string containerPath = appendPrefix(appendPrefix(filteredPath, "properties"), sidreGroup->getName());
  // auto& containerNode = container.name().empty() ? m_schemaRoot : m_schemaRoot[containerPath];

  // m_inletContainerPathNames.push_back(sidreGroup->getPathName());
  // auto& currContainer =
  //   m_rstTables.emplace(sidreGroup->getPathName(), ContainerData {m_colLabels})
  //     .first->second;
  // currContainer.containerName = sidreGroup->getName();
  if(sidreGroup->getName() != "" && sidreGroup->hasView("description"))
  {
    containerNode["description"] =
      std::string(sidreGroup->getView("description")->getString());
  }

  // // FIXME: Handle container fields differently
  for(const auto& fieldEntry : container.getChildFields())
  {
    const auto name = removeBeforeDelimiter(fieldEntry.first);
    auto& childNode = containerNode["properties"][name];
    detail::recordFieldSchema(*fieldEntry.second, childNode);
    if(fieldEntry.second->isRequired())
    {
      containerNode["required"].append() = name;
    }
  }
}

void JSONSchemaWriter::finalize() { m_schemaRoot.save(m_fileName, "json"); }

}  // namespace inlet
}  // namespace axom
