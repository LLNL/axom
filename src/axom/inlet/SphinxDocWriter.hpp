#include <string>
#include <vector>
#include <fstream>

#include "axom/sidre.hpp"
#include "axom/inlet/DocWriter.hpp"

#ifndef INLET_SPHINXDOCWRITER_HPP
#define INLET_SPHINXDOCWRITER_HPP

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
  void writeTable(const std::string& title);

  std::ofstream m_outFile;
  std::ostringstream m_oss;
  axom::sidre::Group* m_sidreRootGroup;
  std::vector<std::vector<std::string>> m_rstTable;
  std::string m_fileName;
};

}
}

#endif