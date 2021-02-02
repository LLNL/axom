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
Table& Inlet::addTable(const std::string& name, const std::string& description)
{
  return m_globalTable.addTable(name, description);
}

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

void Inlet::registerDocWriter(std::unique_ptr<DocWriter> writer)
{
  m_docWriter = std::move(writer);
}

namespace detail
{
/*!
 *******************************************************************************
 * \brief Recursive helper function for traversing an Inlet tree for documentation
 * generation purposes
 * 
 * \param [inout] writer The DocWriter to use for documentation
 * \param [in] table The current table to write
 *******************************************************************************
 */
void docWriterHelper(DocWriter& writer, const Table& table)
{
  // Use a pre-order traversal for readability
  writer.documentTable(table);
  for(const auto& sub_table_entry : table.getChildTables())
  {
    docWriterHelper(writer, *sub_table_entry.second);
  }
}

}  // end namespace detail

void Inlet::writeDoc()
{
  if(m_docEnabled)
  {
    detail::docWriterHelper(*m_docWriter, m_globalTable);
    m_docWriter->finalize();
  }
}

bool Inlet::verify() const { return m_globalTable.verify(); }

}  // end namespace inlet
}  // end namespace axom
