// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file SphinxDocWriter.cpp
 *
 * \brief This file contains the class implementation of the SphinxDocWriter.
 *******************************************************************************
 */

#include "axom/inlet/SphinxDocWriter.hpp"
#include <iostream>
#include <assert.h>
#include "axom/slic.hpp"

namespace axom
{
namespace inlet
{
 
SphinxDocWriter::SphinxDocWriter(const std::string& fileName, axom::sidre::Group* root) {
  SLIC_WARNING_IF(root != nullptr, "Sidre Group is null pointer");
  m_sidreRootGroup = root;
  m_fileName = fileName;
}

void SphinxDocWriter::writeDocumentation() {
  if (m_sidreRootGroup->getName() == "") {
    writeTitle("Input Deck Options");
  }  else {
    writeTitle(m_sidreRootGroup->getName());
  }
  m_oss << ".. |uncheck|    unicode:: U+2610 .. UNCHECKED BOX" << std::endl;
  m_oss << ".. |check|      unicode:: U+2611 .. CHECKED BOX" << std::endl;
  writeDocumentationHelper(m_sidreRootGroup);
  m_rstTables.push_back(m_currentTable);
  writeAllTables();
  m_outFile.open(m_fileName);
  m_outFile << m_oss.str();
  m_outFile.close();
}

void SphinxDocWriter::writeDocumentationHelper(axom::sidre::Group* sidreGroup) {
  SLIC_WARNING_IF(sidreGroup, "Root was nullptr");
  axom::sidre::IndexType i = sidreGroup->getFirstValidGroupIndex();

  // Case 1: the current group is a Field so attributes are stored in views

  if (sidreGroup != m_sidreRootGroup && i == axom::sidre::InvalidIndex) { 
    collectFieldInfo(sidreGroup);
  } 

  // Case 2: Current root corresponds to a inlet::Table

  if (i != axom::sidre::InvalidIndex) {
    m_rstTables.push_back(m_currentTable);
    m_currentTable = TableData();
    m_currentTable.tableName = sidreGroup->getName();
    if (sidreGroup->getName() != "" && sidreGroup->hasChildView("description")) {
      m_currentTable.description = sidreGroup->getView("description")->getString();
    } 
  }

  while (i != axom::sidre::InvalidIndex) {
    writeDocumentationHelper(sidreGroup->getGroup(i));
    i = sidreGroup->getNextValidGroupIndex(i);
  }
}

void SphinxDocWriter::writeTitle(const std::string& title) {
  if (title != "") {
    std::string equals = std::string(title.length(), '=');
    m_oss << equals << "\n" << title << "\n" << equals << "\n";
  }
}

void SphinxDocWriter::writeSubtitle(const std::string& sub) {
  if (sub != "") {
    std::string dashes = std::string(sub.length(), '-');
    m_oss << "\n" << dashes << "\n" << sub << "\n" << dashes << "\n\n";
  }
}

void SphinxDocWriter::writeTable(const std::string& title, 
                                 const std::vector<std::vector<std::string>>& rstTable) {
  SLIC_WARNING_IF(rstTable.size() > 1, "Vector for corresponding rst table must be nonempty");
  std::string result = ".. list-table:: " + title;
  result += "\n   :widths: 25 25 25 25 25\n   :header-rows: 1\n   :stub-columns: 1\n\n";
  for (unsigned int i = 0; i < rstTable.size(); i++) {
    result += "   * - ";
    for (unsigned int j = 0; j < rstTable[i].size(); j++) {
      if (j != 0) {
        result += "     - ";
      }
      result += rstTable[i][j] + "\n";
    }
  }
  m_oss << result;
}

void SphinxDocWriter::writeAllTables() {
  for (unsigned int i = 0; i < m_rstTables.size(); i++) {
    writeSubtitle(m_rstTables[i].tableName);
    if (m_rstTables[i].description != "") {
      m_oss << "Description: " << m_rstTables[i].description << std::endl << std::endl;
    }
    if (m_rstTables[i].rstTable.size() > 1) {
      writeTable("Fields", m_rstTables[i].rstTable);
    }
  }

}

void SphinxDocWriter::collectFieldInfo(axom::sidre::Group* sidreGroup) {
  std::vector<std::string> fieldAttributes(5, "");
  fieldAttributes.resize(5);

  fieldAttributes[0] = sidreGroup->getName();
  if (sidreGroup->hasView("description")) {
    fieldAttributes[1] = std::string(sidreGroup->getView("description")->getString());
  } 

  if (sidreGroup->hasView("defaultValue")) {
    axom::sidre::TypeID type = sidreGroup->getView("defaultValue")->getTypeID();
    if (type == axom::sidre::TypeID::INT8_ID) {
      int8 val = sidreGroup->getView("defaultValue")->getData();
      fieldAttributes[2] = std::to_string(val);
    } else if (type == axom::sidre::TypeID::INT_ID) {
      int val = sidreGroup->getView("defaultValue")->getData();
      fieldAttributes[2] = std::to_string(val);
    } else if (type == axom::sidre::DOUBLE_ID) {
      double val = sidreGroup->getView("defaultValue")->getData();
      fieldAttributes[2] = std::to_string(val);
    } else {
      fieldAttributes[2] = sidreGroup->getView("defaultValue")->getString();
    }
  }

  if (sidreGroup->hasView("range")) {
    axom::sidre::TypeID type = sidreGroup->getView("range")->getTypeID();
    if (type == axom::sidre::INT_ID) {
      int* range = sidreGroup->getView("range")->getArray();
      fieldAttributes[3] = std::to_string(range[0]) + " to " + std::to_string(range[1]);
    } else  if (type == axom::sidre::DOUBLE_ID) {
      double* range = sidreGroup->getView("range")->getArray();
      fieldAttributes[3] = std::to_string(range[0]) + " to " + std::to_string(range[1]);
    } else {
      SLIC_WARNING("Unexpected range type");
    }
  } else if (sidreGroup->hasView("validValues")) {
    SLIC_WARNING_IF(sidreGroup->getView("validValues")->getTypeID() == axom::sidre::INT_ID,
                    "discrete range is only valid for integers");
    int* range = sidreGroup->getView("validValues")->getArray();
    size_t size = sidreGroup->getView("validValues")->getBuffer()->getNumElements();
    for (size_t i = 0; i < size; i++) {
      if (i == size-1) {
        fieldAttributes[3] += std::to_string(range[i]);      
      } else {
        fieldAttributes[3] += std::to_string(range[i]) + ", ";
      }
    }
  }
  
  if (sidreGroup->hasView("required")) {
    int8 required = sidreGroup->getView("required")->getData();
    fieldAttributes[4] = required ? "|check|" : "|uncheck|";
  } else {
    fieldAttributes[4] = "|uncheck|";
  }

  m_currentTable.rstTable.push_back(fieldAttributes);
}

}
}
