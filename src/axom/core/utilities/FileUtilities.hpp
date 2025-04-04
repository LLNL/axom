// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef COMMON_FILE_UTILITIES_H_
#define COMMON_FILE_UTILITIES_H_

#include <string>

namespace axom
{
namespace utilities
{
namespace filesystem
{
/*!
 * \brief Returns the current working directory within the file system
 * \return path of cwd if found, else, returns "./"
 */
std::string getCWD();

/*!
 * \brief Changes the current working directory within the file system
 * \param [in] dirName an absolute or relative directory path
 * \return Status code 0 for success, non-zero otherwise
 */
int changeCWD(const std::string& dirName);

/*!
 * \brief Checks if the path pointed to by fileName exists
 * \param [in] fileName name of the file to check
 * \return true if file system contains a file named fileName, false otherwise
 */
bool pathExists(const std::string& fileName);

/*!
 * \brief Joins two strings with the given separator char
 *
 * \param [in] fileDir The directory of the file
 * \param [in] fileName The name of the file
 * \param [in] separator a single character to seperate the two strings
 *
 * \note fileDir can be the empty string
 * \note fileName can include directories (e.g. a/b/c.txt), but should not be an
 * absolute path f a non-empty fileDir is supplied.
 *
 * \returns The concatenated string: fileDir + fileName, with separator in
 *  between, if necessary
 * \note Example1:  joinPath("abc", "def") -> "abc/def"
 * \note Example2:  joinPath("abc/", "def") -> "abc/def"
 * \note Example3:  joinPath("abc/", "def/ghi") -> "abc/def/ghi"
 */
std::string joinPath(const std::string& fileDir,
                     const std::string& fileName,
                     const std::string& separator = "/");

/*!
 * \brief Make directories for a given path string
 *
 * \param [in] path  string representing an absolute or relative directory path
 *
 * Everything in the path is assumed to be intended to be a directory.  If
 * a directory in the path already exists, nothing is done.  If a directory
 * doesn't exist, it is created.
 */
int makeDirsForPath(const std::string& path);

/*!
 * \brief Add a prefix to a path if it is relative.
 *
 * \param [in] path non-empty string representing an absolute or relative path
 * \param [in] prefix string representing a directory path
 *
 * \return \c prefix + \c path (with delimiter) if the path is
 * relative or \c prefix is empty, or just \c path otherwise.
 */
std::string prefixRelativePath(const std::string& path, const std::string& prefix);

/*!
 * \brief Get parent path name from a filesystem path.
 *
 * \param [in] path a non-empty absolute or relative filesystem path
 *
 * \return a directory path formed by removing the last part of the
 *    input path
 *
 * If path is the system root, its parent is the empty string.
 */
std::string getParentPath(const std::string& path);

/*!
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
 * string will be "abc/def/ghi".  If the path string is "file.txt", the
 * output dir string will be "" (the empty string).
 */
void getDirName(std::string& dir, const std::string& path);

/*!
 * \brief Remove the specified file.
 * \param filename The name of the file.
 * \return 0 on success, -1 on failure. errno can obtain more information
 *         about the failure.
 *
 * \note On Windows, this function calls _unlink() and will fail if there are
 *       any open file handles to the specified file.
 *       On Linux, this function calls unlink() and will succeed even if
 *       there are open file handles to the specified file.
 */
int removeFile(const std::string& filename);

}  // end namespace filesystem
}  // end namespace utilities
}  // end namespace axom

#endif  //  COMMON_FILE_UTILITIES_H_
