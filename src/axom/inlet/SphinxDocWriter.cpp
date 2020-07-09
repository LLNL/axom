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
  m_outFile.open(m_fileName);
  if (sidreGroup->getName() == "") {
    writeTitle("Untitled");
  }  else {
    writeTitle(sidreGroup->getName());
  }
  m_rstTable = {{"Field Name", "Description", "Default Value", "Range", "Required"}};
  writeDocumentsHelper(m_sidreRootGroup);
  writeTable("Fields");
  m_outFile.close();
}

void SphinxDocWriter::writeDocumentsHelper(axom::sidre::Group* sidreGroup) {
  SLIC_ASSERT_MSG(sidreGroup, "Root was nullptr");
  axom::sidre::IndexType i = sidreGroup->getFirstValidGroupIndex();
  if (sidreGroup != m_sidreRootGroup && i == axom::sidre::InvalidIndex) { 
    // means that it is a field so attributes are stored in views
    std::vector<std::string> fieldAttributes(5, "Not specified");
    fieldAttributes.resize(5);
    fieldAttributes[0] = sidreGroup->getName();
    if (sidreGroup->hasView("description")) {
      fieldAttributes[1] = std::string(sidreGroup->getView(sidreGroup->getViewIndex("description"))->getString());
    } 
    // fieldAttributes[2] = std::to_string(sidreGroup->getView(sidreGroup->getViewIndex("default values"))->getScalar());
    // fieldAttributes[3] = std::to_string(sidreGroup->getView(sidreGroup->getViewIndex("range"))->getScalar());
    if (sidreGroup->hasView("required")) {
      int8 required = sidreGroup->getView(sidreGroup->getViewIndex("required"))->getData();
      fieldAttributes[4] = required ? "True" : "False";
    } 
  
    m_rstTable.push_back(fieldAttributes);
  } 

  if (i != axom::sidre::InvalidIndex) {
    writeSubtitle(sidreGroup->getName());
    if (sidreGroup->getName() != "" && sidreGroup->hasChildView("description")) {
      m_outFile << "\nDescription: " << sidreGroup->getView("description")->getString() << "\n\n";
    } 
  }

  while (i != axom::sidre::InvalidIndex) {
    writeDocumentsHelper(sidreGroup->getGroup(i));
    i = sidreGroup->getNextValidGroupIndex(i);
  }
}

void SphinxDocWriter::writeTitle(const std::string& title) {
  SLIC_ASSERT_MSG(m_outFile.is_open(), "Output file should be open");
  if (title != "") {
    std::string equals = std::string(title.length(), '=');
    m_outFile << equals << "\n" << title << "\n" << equals << "\n";
  }
}

void SphinxDocWriter::writeSubtitle(const std::string& sub) {
  SLIC_ASSERT_MSG(m_outFile.is_open(), "Output file should be open");  
  
  if (sub != "") {
    std::string dashes = std::string(sub.length(), '-');
    m_outFile << dashes << "\n" << sub << "\n" << dashes << "\n";
  }
}

void SphinxDocWriter::writeTable(const std::string& title) {
  SLIC_ASSERT_MSG(m_outFile.is_open(), "Output file should be open");  
  std::string result = ".. list-table:: " + title;
  result += "\n   :widths: 25 25 25 25 25\n   :header-rows: 1\n   :stub-columns: 1\n\n";
  for (int i = 0; i < m_rstTable.size(); i++) {
    result += "   * - ";
    for (int j = 0; j < m_rstTable[i].size(); j++) {
      if (j != 0) {
        result += "     - ";
      }
      result += m_rstTable[i][j] + "\n";
    }
  }
  m_outFile << result;
}

}
}