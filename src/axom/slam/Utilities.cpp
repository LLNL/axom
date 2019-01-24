/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

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
