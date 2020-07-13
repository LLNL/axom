// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef INLET_DOCWRITER_HPP
#define INLET_DOCWRITER_HPP

#include <string>
#include <vector>
#include <fstream>

#include "axom/sidre.hpp"

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