#include <string>
#include <vector>
#include <fstream>

#include "axom/sidre.hpp"

#ifndef INLET_DOCWRITER_HPP
#define INLET_DOCWRITER_HPP

namespace axom
{
namespace inlet
{ 

class DocWriter {
public: 
  DocWriter(const std::string& fileName, axom::sidre::Group* sidreRootGroup);
  void writeTitle(const std::string& title);
  void writeSubtitle(const std::string& sub);
  void writeTable(const std::string& title);

  // Entry point for document generation
  void writeDocuments(axom::sidre::Group* sidreGroup);

private:
  std::ofstream outFile;
  axom::sidre::Group* sidreGroupRoot;
  std::vector<std::vector<std::string>> rstTable;
};

}
}

#endif