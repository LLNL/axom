// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file Inlet.cpp
 *
 * \brief This file contains the class implementation of Inlet, the main class
 *        for the Inlet component.
 *******************************************************************************
 */

#include "axom/inlet/Inlet.hpp"

#include "fmt/fmt.hpp"
#include "axom/slic.hpp"
#include "axom/inlet/inlet_utils.hpp"
#include <algorithm>

namespace axom
{
namespace inlet
{
Container& Inlet::addStruct(const std::string& name,
                            const std::string& description)
{
  return m_globalContainer.addStruct(name, description);
}

VerifiableScalar& Inlet::addBool(const std::string& name,
                                 const std::string& description)
{
  return m_globalContainer.addBool(name, description);
}

VerifiableScalar& Inlet::addDouble(const std::string& name,
                                   const std::string& description)
{
  return m_globalContainer.addDouble(name, description);
}

VerifiableScalar& Inlet::addInt(const std::string& name,
                                const std::string& description)
{
  return m_globalContainer.addInt(name, description);
}

VerifiableScalar& Inlet::addString(const std::string& name,
                                   const std::string& description)
{
  return m_globalContainer.addString(name, description);
}

void Inlet::registerWriter(std::unique_ptr<Writer> writer)
{
  m_writer = std::move(writer);
}

namespace detail
{
/*!
 *******************************************************************************
 * \brief Recursive helper function for traversing an Inlet tree for documentation
 * generation purposes
 * 
 * \param [inout] writer The Writer to use for documentation
 * \param [in] container The current container to write
 *******************************************************************************
 */
void writerHelper(Writer& writer, const Container& container)
{
  // Use a pre-order traversal for readability
  writer.documentContainer(container);
  // If the current container contains a collection, visit that last
  const auto& child_containers = container.getChildContainers();
  for(const auto& sub_container_entry : child_containers)
  {
    // Ignore the collection group as it will be visited later
    if(!isCollectionGroup(sub_container_entry.first))
    {
      writerHelper(writer, *sub_container_entry.second);
    }
  }
  auto iter = child_containers.find(
    appendPrefix(container.name(), detail::COLLECTION_GROUP_NAME));
  if(iter != child_containers.end())
  {
    const auto& coll_container = *iter->second;
    if(coll_container.sidreGroup()->hasGroup(detail::COLLECTION_INDICES_NAME))
    {
      writerHelper(writer, coll_container);
    }
  }
}

}  // end namespace detail

void Inlet::writeDoc()
{
  if(m_docEnabled)
  {
    detail::writerHelper(*m_writer, m_globalContainer);
    m_writer->finalize();
  }
}

bool Inlet::verify() const { return m_globalContainer.verify(); }

}  // end namespace inlet
}  // end namespace axom
