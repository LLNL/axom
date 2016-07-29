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

namespace asctoolkit {
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
  std::string joinPath(const std::string& fileDir
                       , const std::string& fileName
                       , const std::string& separator = "/");


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
 * \brief Make directories for a given path string that contains a file name
 *
 * \param [in] path  string representing an absolute or relative directory path
 * \param [in] mode  controls the permissions of the new directories
 *
 * This is similar to makeDirsForPath, exept it assumes that the last word in
 * the path string is intended to be a file name and thus will not create a
 * directory with that name.  For example, if the path string is
 * "abc/def/ghi/file.txt", this function will create the directory path
 * abc/def/ghi (if it doesn't already exist) and ignore the file name at the
 * end of the string.
 */
  int truncatedMakeDirs(
    const std::string& path,
    mode_t mode = 0700);

} // end namespace filesystem
} // end namespace utilities
} // end namespace asctoolkit

#endif //  COMMON_FILE_UTILITIES_H_
