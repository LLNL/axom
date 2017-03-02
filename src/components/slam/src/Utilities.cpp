/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

#include "slam/Utilities.hpp"

#include "common/FileUtilities.hpp"
#include "slic/slic.hpp"

namespace axom {
namespace slam {
namespace util {


  std::string findFileInAncestorDirs(const std::string& fileName, int numAncestors)
  {
    using namespace axom::utilities::filesystem;

    bool fileFound = pathExists( fileName);
    std::string ancestorFileName = fileName;

    std::string ancestorPath = "..";

    // File not found -- try to find path by walking up tree
    for(int i = 0; !fileFound && i< numAncestors; ++i)
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
      SLIC_WARNING( "Could not find file: '"  << fileName
                                              << "' (also tried several ancestor directories')"
                                              << "\nThe current working directory is: '" << getCWD() << "'");

      ancestorFileName = fileName;
    }

    return ancestorFileName;
  }



} // end namespace util
} // end namespace slam
} // end namespace axom
