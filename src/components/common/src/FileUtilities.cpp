/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

/*!
 ******************************************************************************
 *
 * \file
 *
 * \brief   Implementation file for file utility functions.
 *
 ******************************************************************************
 */

#include "common/FileUtilities.hpp"

#include <string>
#include <fstream>
#include <sstream>


#include <cstdio>                       // defines FILENAME_MAX

#ifdef WINDOWS
    #include <direct.h>
    #include <sys/stat.h>
    // Warning: not yet tested on windows
    #define GetCurrentDir _getcwd
    #define Stat _stat
#else
    #include <unistd.h>                 // for getcwd
    #include <sys/stat.h>               // for stat
    #define GetCurrentDir getcwd
    #define Stat stat
#endif



namespace asctoolkit {
namespace utilities {
namespace filesystem {

  std::string getCWD()
  {
    char cCurrentPath[FILENAME_MAX];

    if (!GetCurrentDir(cCurrentPath, FILENAME_MAX))
    {
      //Note: Cannot use logging in COMMON component -- topic of JIRA issue ATK-463
      //SLIC_WARNING("Common::Could not find cwd");

      return std::string("./");
    }

    return std::string(cCurrentPath);
  }

//-----------------------------------------------------------------------------
  bool pathExists(const std::string& fileName)
  {
    // Uses system's stat() function.
    // Return code 0 indicates file exists
    struct Stat buffer;
    return (Stat(fileName.c_str(), &buffer) == 0);
  }

//-----------------------------------------------------------------------------
  std::string joinPath(const std::string& fileDir, const std::string& fileName,
                       const std::string& separator)
  {
    // Check if we need to add a separator
    bool pathNeedsSep = !fileDir.empty()
                        && (fileDir[ fileDir.size() -1] != separator[0]);

    // Concatenate the path with the fileName to create the full path
    std::stringstream fullFileNameStream;
    fullFileNameStream << fileDir
                       << (pathNeedsSep ? separator : "" )
                       << fileName;

    return fullFileNameStream.str();
  }

//-----------------------------------------------------------------------------
  int makeDirsForPath(const std::string& path, mode_t mode)
  {

    char separator = '/';
    std::string::size_type pos = 0;
    int err = 0;

    do {
      pos = path.find(separator, pos+1);
      std::string dir_name = path.substr(0, pos);
      err = mkdir(dir_name.c_str(), mode);
      err = (err && (errno != EEXIST)) ? 1 : 0;

    } while (pos != std::string::npos);

    return err;
  }

//-----------------------------------------------------------------------------
  int truncatedMakeDirs(const std::string& path, mode_t mode)
  {
    int err = 0;

    std::string curr;
    std::string dir_path;
    std::string slash = "/";

    std::size_t found = path.rfind(slash);
    if (found != std::string::npos)
    {
      dir_path = path.substr(0,found);
      if (found != path.size()-1) {
        curr = path.substr(found+1, path.size()-(found-1));
      }
    }
    else
    {
      curr = path;
    }
    if (!dir_path.empty()) {
      err = makeDirsForPath(dir_path, mode);
    }

    return err;
  }

}   // end namespace filesystem
}   // end namespace utilities
}   // end namespace asctoolkit
