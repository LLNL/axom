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

void filterCollectionPaths(std::string& target,
                           const std::vector<std::string>& collectionPaths,
                           const std::string& label)
{
  for(const auto& path : collectionPaths)
  {
    // Needs to be at the beginning - since the path is absolute
    if(target.find(path) == 0)
    {
      // Remove the name of the selected element - get the slash following the name
      const auto pos = target.find('/', path.length() + 1);
      if(pos != std::string::npos)
      {
        target.erase(path.length() + 1, pos - path.length());
      }
      // Insert the extension just after it
      target.insert(path.length(), "/" + label);
    }
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

  // The "level" of the schema to add to - adjusted "up" if we're an element
  // of a collection
  const sidre::Group* propertyGroup = sidreGroup;
  if(isCollectionGroup(sidreGroup->getParent()->getName()))
  {
    propertyGroup = sidreGroup->getParent();
  }
  // FIXME: Use path class for some of this logic - all this is doing is inserting
  // a 'properties' in after each token
  std::string filteredPathName =
    removePrefix(m_rootPath, propertyGroup->getPathName());
  detail::filterCollectionPaths(filteredPathName,
                                m_ArrayPaths,
                                "additionalItems");
  detail::filterCollectionPaths(filteredPathName,
                                m_DictionaryPaths,
                                "additionalProperties");
  std::vector<std::string> tokens;
  utilities::string::split(tokens, filteredPathName, '/');
  // Remove the collection group annotations
  auto iter =
    std::remove(tokens.begin(), tokens.end(), detail::COLLECTION_GROUP_NAME);
  tokens.erase(iter, tokens.end());
  std::string containerPath =
    fmt::format("properties/{}", fmt::join(tokens, "/properties/"));
  // Needed so we can reassign without calling the op= on the Node itself
  std::reference_wrapper<conduit::Node> containerNode =
    container.name().empty() ? m_schemaRoot : m_schemaRoot[containerPath];

  // Annotate with the default type of object, but arrays should be labeled
  // with "array"
  if(!containerNode.get().has_child("type"))
  {
    containerNode.get()["type"] = "object";
  }
  if(isCollectionGroup(container.name()))
  {
    auto indices = detail::collectionIndices(container);
    // Check to see if it's an array - have to watch out for mixed index types
    if(std::all_of(indices.begin(), indices.end(), [](const VariantKey& key) {
         return key.type() == InletType::Integer;
       }))
    {
      containerNode.get()["type"] = "array";
      m_ArrayPaths.push_back(filteredPathName);
    }
    else
    {
      m_DictionaryPaths.push_back(filteredPathName);
    }
  }

  if(sidreGroup->getName() != "" && sidreGroup->hasView("description"))
  {
    containerNode.get()["description"] =
      std::string(sidreGroup->getView("description")->getString());
  }

  for(const auto& fieldEntry : container.getChildFields())
  {
    const auto name = removeBeforeDelimiter(fieldEntry.first);
    auto& childNode = containerNode.get()["properties"][name];
    detail::recordFieldSchema(*fieldEntry.second, childNode);
    if(fieldEntry.second->isRequired())
    {
      containerNode.get()["required"].append() = name;
    }
  }
}

void JSONSchemaWriter::finalize() { m_schemaRoot.save(m_fileName, "json"); }

}  // namespace inlet
}  // namespace axom
