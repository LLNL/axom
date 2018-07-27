/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
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

#include "axom/core/utilities/FileUtilities.hpp"

#include <string>
#include <fstream>
#include <sstream>
#include <cerrno>

#include <cstdio>                       // defines FILENAME_MAX

#ifdef WIN32
    #include <direct.h>
    #include <sys/stat.h>

    #define GetCurrentDir _getcwd
    #define Stat _stat
#else
    #include <unistd.h>                 // for getcwd
    #include <sys/stat.h>               // for stat
    #define GetCurrentDir getcwd
    #define Stat stat
#endif



namespace axom
{
namespace utilities
{
namespace filesystem
{

std::string getCWD()
{
  char cCurrentPath[FILENAME_MAX];

  if (!GetCurrentDir(cCurrentPath, FILENAME_MAX))
  {
    //Note: Cannot use logging in COMMON component -- topic of JIRA issue
    // ATK-463
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
int makeDirsForPath(const std::string& path)
{

  char separator = '/';
  std::string::size_type pos = 0;
  int err = 0;

  do
  {
    pos = path.find(separator, pos+1);
    std::string dir_name = path.substr(0, pos);
#ifdef WIN32
    err = _mkdir(dir_name.c_str());
#else
    mode_t mode = 0770;     // user and group rwx permissions
    err = mkdir(dir_name.c_str(), mode);
#endif
    err = (err && (errno != EEXIST)) ? 1 : 0;

  }
  while (pos != std::string::npos);

  return err;
}

//-----------------------------------------------------------------------------
void getDirName(std::string& dir, const std::string& path)
{
  char separator = '/';

  std::size_t found = path.rfind(separator);
  if (found != std::string::npos)
  {
    dir = path.substr(0,found);
  }
  else
  {
    dir = "";
  }
}

}   // end namespace filesystem
}   // end namespace utilities
}   // end namespace axom
