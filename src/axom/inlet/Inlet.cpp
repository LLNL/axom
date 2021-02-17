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
Table& Inlet::addStruct(const std::string& name, const std::string& description)
{
  return m_globalTable.addStruct(name, description);
}

VerifiableScalar& Inlet::addBool(const std::string& name,
                                 const std::string& description)
{
  return m_globalTable.addBool(name, description);
}

VerifiableScalar& Inlet::addDouble(const std::string& name,
                                   const std::string& description)
{
  return m_globalTable.addDouble(name, description);
}

VerifiableScalar& Inlet::addInt(const std::string& name,
                                const std::string& description)
{
  return m_globalTable.addInt(name, description);
}

VerifiableScalar& Inlet::addString(const std::string& name,
                                   const std::string& description)
{
  return m_globalTable.addString(name, description);
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
 * \param [in] table The current table to write
 *******************************************************************************
 */
void writerHelper(Writer& writer, const Table& table)
{
  // Use a pre-order traversal for readability
  writer.documentTable(table);
  // If the current table contains a collection, visit that last
  const auto& child_tables = table.getChildTables();
  for(const auto& sub_table_entry : child_tables)
  {
    // Ignore the collection group as it will be visited later
    if(!isCollectionGroup(sub_table_entry.first))
    {
      writerHelper(writer, *sub_table_entry.second);
    }
  }
  auto iter =
    child_tables.find(appendPrefix(table.name(), detail::COLLECTION_GROUP_NAME));
  if(iter != child_tables.end())
  {
    const auto& coll_table = *iter->second;
    if(coll_table.sidreGroup()->hasGroup(detail::COLLECTION_INDICES_NAME))
    {
      writerHelper(writer, coll_table);
    }
  }
}

}  // end namespace detail

void Inlet::writeDoc()
{
  if(m_docEnabled)
  {
    detail::writerHelper(*m_writer, m_globalTable);
    m_writer->finalize();
  }
}

bool Inlet::verify() const { return m_globalTable.verify(); }

}  // end namespace inlet
}  // end namespace axom
