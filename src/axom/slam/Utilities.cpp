// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/slam/Utilities.hpp"

#include "axom/core/utilities/FileUtilities.hpp"
#include "axom/slic/interface/slic.hpp"

namespace axom
{
namespace slam
{
namespace util
{


std::string
findFileInAncestorDirs(const std::string& fileName, int numAncestors)
{
  using namespace axom::utilities::filesystem;

  bool fileFound = pathExists( fileName);
  std::string ancestorFileName = fileName;

  std::string ancestorPath = "..";

  // File not found -- try to find path by walking up tree
  for(int i = 0 ; !fileFound && i< numAncestors ; ++i)
  {
    ancestorFileName = joinPath(ancestorPath, fileName);
    fileFound = pathExists( ancestorFileName );

    if( !fileFound )
    {
      ancestorPath = "../" + ancestorPath;
    }
  }

  if(!fileFound)
  {
    SLIC_WARNING(
      "Could not find file: '"
      << fileName << "' (also tried several ancestor directories')"
      << "\nThe current working directory is: '" << getCWD() << "'");

    ancestorFileName = fileName;
  }

  return ancestorFileName;
}



} // end namespace util
} // end namespace slam
} // end namespace axom
