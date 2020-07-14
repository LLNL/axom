// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef INLET_SPHINXDOCWRITER_HPP
#define INLET_SPHINXDOCWRITER_HPP

#include <string>
#include <vector>
#include <fstream>

#include "axom/sidre.hpp"
#include "axom/inlet/DocWriter.hpp"

namespace axom
{
namespace inlet
{ 

class SphinxDocWriter : public DocWriter {

public: 
  SphinxDocWriter(const std::string& fileName, axom::sidre::Group* sidreRootGroup);

  // Entry point for document generation
  void writeDocuments(axom::sidre::Group* sidreGroup);

private:
  void writeDocumentsHelper(axom::sidre::Group* sidreGroup);

  void writeTitle(const std::string& title);
  void writeSubtitle(const std::string& sub);
  void writeTable(const std::string& title, const std::vector<std::vector<std::string>>& rstTable);
  void writeAllTables();

  struct TableData {
    TableData() {
      rstTable = {{"Field Name", "Description", "Default Value", "Range", "Required"}};
    }
    std::string tableName;
    std::string description;
    std::vector<std::vector<std::string>> rstTable;
  };

  std::ofstream m_outFile;
  std::ostringstream m_oss;
  axom::sidre::Group* m_sidreRootGroup;
  std::vector<TableData> m_rstTables;
  std::string m_fileName;
  TableData m_currentTable;
};

}
}

#endif