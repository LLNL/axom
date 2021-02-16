// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file SphinxWriter.cpp
 *
 * \brief This file contains the class implementation of the SphinxWriter.
 *******************************************************************************
 */

#include "axom/inlet/SphinxWriter.hpp"

#include <iostream>

#include "axom/slic.hpp"
#include "axom/inlet/Table.hpp"

namespace axom
{
namespace inlet
{
SphinxWriter::SphinxWriter(const std::string& fileName)
  : m_colLabels({"Field Name",
                 "Description",
                 "Default Value",
                 "Range/Valid Values",
                 "Required"})
{
  m_fileName = fileName;

  m_oss << ".. |uncheck|    unicode:: U+2610 .. UNCHECKED BOX\n";
  m_oss << ".. |check|      unicode:: U+2611 .. CHECKED BOX\n\n";
  writeTitle("Input file Options");
}

void SphinxWriter::documentTable(const Table& table)
{
  const auto sidreGroup = table.sidreGroup();
  m_inletTablePathNames.push_back(sidreGroup->getPathName());
  auto& currTable =
    m_rstTables.emplace(sidreGroup->getPathName(), TableData {m_colLabels})
      .first->second;
  currTable.tableName = sidreGroup->getName();
  if(sidreGroup->getName() != "" && sidreGroup->hasView("description"))
  {
    currTable.description = sidreGroup->getView("description")->getString();
  }

  // FIXME: Handle container fields differently
  for(const auto& field_entry : table.getChildFields())
  {
    extractFieldMetadata(field_entry.second->sidreGroup());
  }
}

void SphinxWriter::finalize()
{
  writeAllTables();
  m_outFile.open(m_fileName);
  m_outFile << m_oss.str();
  m_outFile.close();
}

void SphinxWriter::writeTitle(const std::string& title)
{
  if(title != "")
  {
    std::string equals = std::string(title.length(), '=');
    m_oss << equals << "\n" << title << "\n" << equals << "\n";
  }
}

void SphinxWriter::writeSubtitle(const std::string& sub)
{
  if(sub != "")
  {
    std::string dashes = std::string(sub.length(), '-');
    m_oss << "\n" << dashes << "\n" << sub << "\n" << dashes << "\n\n";
  }
}

void SphinxWriter::writeTable(const std::string& title,
                              const std::vector<std::vector<std::string>>& rstTable)
{
  SLIC_WARNING_IF(
    rstTable.size() <= 1,
    "[Inlet] Vector for corresponding rst table must be nonempty");
  std::string result = ".. list-table:: " + title;
  std::string widths = ":widths:";
  // This would be easier with an iterator adaptor like back_inserter but for
  // concatenation
  for(std::size_t i = 0u; i < m_colLabels.size(); i++)
  {
    widths += " 25";
  }
  result += "\n   " + widths + "\n";
  result += "   :header-rows: 1\n   :stub-columns: 1\n\n";
  for(unsigned int i = 0; i < rstTable.size(); ++i)
  {
    result += "   * - ";
    for(unsigned int j = 0; j < rstTable[i].size(); ++j)
    {
      if(j != 0)
      {
        result += "     - ";
      }
      result += rstTable[i][j] + "\n";
    }
  }
  m_oss << result;
}

void SphinxWriter::writeAllTables()
{
  for(std::string& pathName : m_inletTablePathNames)
  {
    auto& currTable = m_rstTables.at(pathName);
    writeSubtitle(currTable.tableName);
    if(currTable.description != "")
    {
      m_oss << "Description: " << currTable.description << std::endl
            << std::endl;
    }
    if(currTable.rstTable.size() > 1)
    {
      writeTable("Fields", currTable.rstTable);
    }
  }
}

std::string SphinxWriter::getValueAsString(const axom::sidre::View* view)
{
  axom::sidre::TypeID type = view->getTypeID();
  if(type == axom::sidre::TypeID::INT8_ID)
  {
    int8 val = view->getData();
    return val ? "True" : "False";
  }
  else if(type == axom::sidre::TypeID::INT_ID)
  {
    int val = view->getData();
    return std::to_string(val);
  }
  else if(type == axom::sidre::TypeID::DOUBLE_ID)
  {
    double val = view->getData();
    return std::to_string(val);
  }
  return view->getString();
}

std::string SphinxWriter::getRangeAsString(const axom::sidre::View* view)
{
  std::ostringstream oss;
  oss.precision(3);
  oss << std::scientific;

  axom::sidre::TypeID type = view->getTypeID();
  if(type == axom::sidre::INT_ID)
  {
    const int* range = view->getData();
    oss << range[0] << " to " << range[1];
  }
  else
  {
    const double* range = view->getData();
    oss << range[0] << " to " << range[1];
  }
  return oss.str();
}

std::string SphinxWriter::getValidValuesAsString(const axom::sidre::View* view)
{
  const int* range = view->getData();
  size_t size = view->getBuffer()->getNumElements();
  std::string result = "";
  for(size_t i = 0; i < size; ++i)
  {
    if(i == size - 1)
    {
      result += std::to_string(range[i]);
    }
    else
    {
      result += std::to_string(range[i]) + ", ";
    }
  }
  return result;
}

std::string SphinxWriter::getValidStringValues(const axom::sidre::Group* sidreGroup)
{
  auto idx = sidreGroup->getFirstValidViewIndex();
  std::string validValues = "";
  while(axom::sidre::indexIsValid(idx))
  {
    validValues += std::string(sidreGroup->getView(idx)->getString());
    idx = sidreGroup->getNextValidViewIndex(idx);
    if(axom::sidre::indexIsValid(idx))
    {
      validValues += ", ";
    }
  }
  return validValues;
}

void SphinxWriter::extractFieldMetadata(const axom::sidre::Group* sidreGroup)
{
  TableData& currentTable =
    m_rstTables.at(sidreGroup->getParent()->getPathName());
  std::vector<std::string> fieldAttributes(m_colLabels.size());

  fieldAttributes[0] = sidreGroup->getName();

  if(sidreGroup->hasView("description"))
  {
    fieldAttributes[1] =
      std::string(sidreGroup->getView("description")->getString());
  }

  if(sidreGroup->hasView("defaultValue"))
  {
    fieldAttributes[2] = getValueAsString(sidreGroup->getView("defaultValue"));
  }

  if(sidreGroup->hasView("range"))
  {
    fieldAttributes[3] = getRangeAsString(sidreGroup->getView("range"));
  }
  else if(sidreGroup->hasView("validValues"))
  {
    fieldAttributes[3] =
      getValidValuesAsString(sidreGroup->getView("validValues"));
  }
  else if(sidreGroup->hasGroup("validStringValues"))
  {
    fieldAttributes[3] =
      getValidStringValues(sidreGroup->getGroup("validStringValues"));
  }

  if(sidreGroup->hasView("required"))
  {
    int8 required = sidreGroup->getView("required")->getData();
    fieldAttributes[4] = required ? "|check|" : "|uncheck|";
  }
  else
  {
    fieldAttributes[4] = "|uncheck|";
  }

  currentTable.rstTable.push_back(fieldAttributes);
}

}  // namespace inlet
}  // namespace axom
