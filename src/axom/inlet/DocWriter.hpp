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
// Entry point for document generation
  virtual void writeDocuments(axom::sidre::Group* sidreGroup) = 0;

private:
  axom::sidre::Group* m_sidreRootGroup;

};

}
}

#endif