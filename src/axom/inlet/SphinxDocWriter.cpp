// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet/SphinxDocWriter.hpp"
#include <iostream>
#include <assert.h>
#include "axom/slic.hpp"

namespace axom
{
namespace inlet
{
 
SphinxDocWriter::SphinxDocWriter(const std::string& fileName, axom::sidre::Group* root) {
  SLIC_ASSERT_MSG(root != nullptr, "Sidre Group is null pointer");
  m_sidreRootGroup = root;
  m_fileName = fileName;
}

void SphinxDocWriter::writeDocuments(axom::sidre::Group* sidreGroup) {
  if (sidreGroup->getName() == "") {
    writeTitle("Input Deck Options");
  }  else {
    writeTitle(sidreGroup->getName());
  }
  m_oss << ".. |uncheck|    unicode:: U+2610 .. UNCHECKED BOX" << std::endl;
  m_oss << ".. |check|      unicode:: U+2611 .. CHECKED BOX" << std::endl;
  writeDocumentsHelper(m_sidreRootGroup);
  m_rstTables.push_back(m_currentTable);
  writeAllTables();
  m_outFile.open(m_fileName);
  m_outFile << m_oss.str();
  m_outFile.close();
}

void SphinxDocWriter::writeDocumentsHelper(axom::sidre::Group* sidreGroup) {
  SLIC_ASSERT_MSG(sidreGroup, "Root was nullptr");
  axom::sidre::IndexType i = sidreGroup->getFirstValidGroupIndex();

  // Case 1: the current group is a Field so attributes are stored in views

  if (sidreGroup != m_sidreRootGroup && i == axom::sidre::InvalidIndex) { 
    std::vector<std::string> fieldAttributes(5, "");
    fieldAttributes.resize(5);
    
    fieldAttributes[0] = sidreGroup->getName();
    if (sidreGroup->hasView("description")) {
      fieldAttributes[1] = std::string(sidreGroup->getView(sidreGroup->getViewIndex("description"))->getString());
    } 

    // the following 2 lines will be needed once default values and range are added to inlet fields
    // fieldAttributes[2] = std::to_string(sidreGroup->getView(sidreGroup->getViewIndex("default values"))->getScalar());
    // fieldAttributes[3] = std::to_string(sidreGroup->getView(sidreGroup->getViewIndex("range"))->getScalar());

    if (sidreGroup->hasView("required")) {
      int8 required = sidreGroup->getView(sidreGroup->getViewIndex("required"))->getData();
      fieldAttributes[4] = required ? "|check|" : "|uncheck|";
    } else {
      fieldAttributes[4] = "|uncheck|";
    }
  
    m_currentTable.rstTable.push_back(fieldAttributes);
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
    writeDocumentsHelper(sidreGroup->getGroup(i));
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
  SLIC_ASSERT_MSG(rstTable.size() > 1, "Vector for corresponding rst table must be nonempty");
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

}
}