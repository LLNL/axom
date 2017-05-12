/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


#ifndef COMMON_FILE_UTILITIES_H_
#define COMMON_FILE_UTILITIES_H_

#include <string>
#include <sys/stat.h>

namespace axom {
namespace utilities {
namespace filesystem {

/**
 * \brief Helper function to print out the current working directory within the file system
 * \return path of cwd if found, else, returns "./"
 */
  std::string getCWD();


/**
 * \brief Helper function to check if the path pointed to by fileName exists
 * \param [in] fileName string name of a file (possibly including relative or absolute path)
 * \return true if file system contains a file named fileName, false otherwise
 */
  bool pathExists(const std::string& fileName);

/**
 * \brief Joins a file directory fileDir and a file name fileName with the given separator char
 * \param [in] fileDir The directory of the file
 * \param [in] fileName The name of the file
 * \note fileDir can be the empty string
 * \note fileName can include directories (e.g. a/b/c.txt), but should not be an absolute path
 *       if a non-empty fileDir is supplied.
 * \param [in] separator The path separator for the file system.  Should be a single character
 * \returns The concatenated string: fileDir + fileName, with separator in between, if necessary
 * \note Example1:  joinPath("abc", "def") -> "abc/def"
 * \note Example2:  joinPath("abc/", "def") -> "abc/def"
 * \note Example3:  joinPath("abc/", "def/ghi") -> "abc/def/ghi"
 */
  std::string joinPath(const std::string& fileDir,
                       const std::string& fileName,
                       const std::string& separator = "/");


/**
 * \brief Make directories for a given path string
 *
 * \param [in] path  string representing an absolute or relative directory path
 * \param [in] mode  controls the permissions of the new directories
 *
 * Everything in the path is assumed to be intended to be a directory.  If
 * a directory in the path already exists, nothing is done.  If a directory
 * doesn't exist, it is created.
 */
  int makeDirsForPath(
    const std::string& path,
    mode_t mode = 0700);

/**
 * \brief Get directory name from a path that contains a file name
 *
 * \param [out] dir  a directory path formed by removing the file name from
 *                   the input path
 * \param [in] path  an absolute or relative directory/file path
 *
 * This function assumes that the input path has a file name at the end, and
 * it removes that file name, leaving a string containing only a directory
 * path.
 * 
 * For example, if the path string is "abc/def/ghi/file.txt", the output dir
 * string will be "abc/def/ghi".
 */
  void getDirName(
    std::string& dir,
    const std::string& path);

} // end namespace filesystem
} // end namespace utilities
} // end namespace axom

#endif //  COMMON_FILE_UTILITIES_H_
