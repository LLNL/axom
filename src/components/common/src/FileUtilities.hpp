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

namespace asctoolkit
{
namespace utilities
{
namespace filesystem
{

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



} // end namespace filesystem
} // end namespace utilities
} // end namespace asctoolkit

#endif //  COMMON_FILE_UTILITIES_H_
