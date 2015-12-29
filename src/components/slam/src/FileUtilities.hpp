/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


#ifndef SLAM_FILE_UTILITIES_H_
#define SLAM_FILE_UTILITIES_H_

#include <string>
#include <fstream>

#include <cstdio>  /* defines FILENAME_MAX */
#ifdef WINDOWS
    #include <direct.h>
    #define GetCurrentDir _getcwd
#else
    #include <unistd.h>
    #define GetCurrentDir getcwd
#endif

#include "slic/slic.hpp"


namespace asctoolkit {
namespace slam {
namespace util {

/**
 * \brief Helper function to find the current working directory within the file system
 * \return a std::string containing the cwd path or a warning message if it cannot be found.
 * \note Uses a SLIC_WARNING to log problems with finding the cwd
 */
  std::string getCWD()
  {
    char cCurrentPath[FILENAME_MAX];

    bool cwdWorked = GetCurrentDir(cCurrentPath, FILENAME_MAX);

    if( !cwdWorked)
    {
        const std::string warnMsg = "< cwd command did not work! >";
        SLIC_WARNING(warnMsg);
        return warnMsg;
    }

    return std::string(cCurrentPath);
  }

  /**
   * \brief Helper function to find a file.
   * If the file does not exist, try to walk up the file system numAncestorsAttempts times.
   * Throws a SLIC error if we cannot find the file after the specified number of times.
   */
  std::string findFileRecursive( const std::string baseFileName
                     , const std::string filePath = ""
                     , int numAncestorsAttempts = 3)
  {
      SLIC_ASSERT(numAncestorsAttempts > 0 );

      std::string origFileName = filePath + baseFileName;

      std::string fileName = origFileName;
      std::ifstream fileStream( fileName.c_str() );

      for(int attempts = 0
              ; !fileStream && attempts < numAncestorsAttempts
              ; ++attempts)
      {
        fileName = "../" + fileName;
        fileStream.open( fileName.c_str());
      }


      SLIC_ERROR_IF( !fileStream
          , "fstream error -- problem opening file: '"  << origFileName
                                                        << "' (also tried several"
                                                        << " ancestors up to '../" << fileName << "')."
                                                        << "\nThe current working directory is: '"
                                                        << getCWD() << "'");
      fileStream.close();

      return fileName;

  }



} // end namespace util
} // end namespace slam
} // end namespace asctoolkit

#endif //  SLAM_FILE_UTILITIES_H_
