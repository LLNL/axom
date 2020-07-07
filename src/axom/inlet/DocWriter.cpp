#include "axom/inlet/DocWriter.hpp"
#include <iostream>
#include <assert.h>
#include "axom/slic.hpp"

namespace axom
{
namespace inlet
{
 
DocWriter::DocWriter(const std::string& fileName, axom::sidre::Group* root, bool isVerbose) {
  SLIC_ASSERT_MSG(root != nullptr, "Sidre Group is null pointer");
  outFile.open(fileName);
  sidreGroupRoot = root;
  verbose = isVerbose;

  if (root->getName() == "") {
    writeTitle("Untitled");
  }  else {
    writeTitle(root->getName());
  }
  rstTable = {{"Field Name", "Description", "Default Value", "Range", "Required"}};
  writeDocuments(sidreGroupRoot);
  writeTable("Fields");
  outFile.close();
}

void DocWriter::writeDocuments(axom::sidre::Group* sidreGroup) {
  SLIC_ASSERT_MSG(sidreGroup, "Root passed into writeDocuments shouldn't be nullptr");
  axom::sidre::IndexType i = sidreGroup->getFirstValidGroupIndex();
  if (sidreGroup != sidreGroupRoot && i == axom::sidre::InvalidIndex) { 
    // means that it is a field so attributes are stored in views
    std::vector<std::string> fieldAttributes(5, "Not specified");
    fieldAttributes.resize(5);
    fieldAttributes[0] = sidreGroup->getName();
    if (verbose) {
      std::cout << "For Field named " << fieldAttributes[0] << ": ";
    }
    if (sidreGroup->hasView("description")) {
      fieldAttributes[1] = std::string(sidreGroup->getView(sidreGroup->getViewIndex("description"))->getString());
      if (verbose) {
        std::cout << "View description found: " << fieldAttributes[1] << ".";
      }
    } else if (verbose) {
      std::cout << "View description not found.";
    }
    // fieldAttributes[2] = std::to_string(sidreGroup->getView(sidreGroup->getViewIndex("default values"))->getScalar());
    // fieldAttributes[3] = std::to_string(sidreGroup->getView(sidreGroup->getViewIndex("range"))->getScalar());
    if (sidreGroup->hasView("required")) {
      int8 required = sidreGroup->getView(sidreGroup->getViewIndex("required"))->getData();
      fieldAttributes[4] = required ? "True" : "False";
      if (verbose) {
        std::cout << " View required found: " << fieldAttributes[4] << ".\n";
      }
    } else if (verbose) {
      std::cout << " View required not found" << std::endl;
    }
  
    rstTable.push_back(fieldAttributes);
  } 

  if (i != axom::sidre::InvalidIndex) {
    writeSubtitle(sidreGroup->getName());
    if (sidreGroup->getName() != "") {
      outFile << "\nDescription: " << sidreGroup->getView("description")->getString() << "\n\n";
      if (verbose) {
        std::cout << "For table named " << sidreGroup->getName() << std::endl;
      }
    } else if (verbose) {
      std::cout << "For base table " << std::endl;
    }
  }

  while (i != axom::sidre::InvalidIndex) {
    writeDocuments(sidreGroup->getGroup(i));
    i = sidreGroup->getNextValidGroupIndex(i);
  }
}

void DocWriter::writeTitle(const std::string& title) {
  SLIC_ASSERT_MSG(outFile.is_open(), "Output file should be open");
  if (title != "") {
    std::string equals;
    for (int i = 0; i < title.length(); i++) {
      equals += "=";
    }
    outFile << equals << "\n" << title << "\n" << equals << "\n";
  }
}

void DocWriter::writeSubtitle(const std::string& sub) {
  SLIC_ASSERT_MSG(outFile.is_open(), "Output file should be open");  
  std::string dashes;
  if (sub != "") {
    for (int i = 0; i < sub.length(); i++) {
      dashes += "-";
    }
    outFile << dashes << "\n" << sub << "\n" << dashes << "\n";
  }
}

void DocWriter::writeTable(const std::string& title) {
  SLIC_ASSERT_MSG(outFile.is_open(), "Output file should be open");  
  std::string result = ".. list-table:: " + title;
  result += "\n   :widths: 25 25 25 25 25\n   :header-rows: 1\n   :stub-columns: 1\n\n";
  for (int i = 0; i < rstTable.size(); i++) {
    result += "   * - ";
    for (int j = 0; j < rstTable[i].size(); j++) {
      if (j != 0) {
        result += "     - ";
      }
      result += rstTable[i][j] + "\n";
    }
  }
  outFile << result;
}

}
}