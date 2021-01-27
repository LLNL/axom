// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef CORE_PATH_UTILITIES_H_
#define CORE_PATH_UTILITIES_H_

#include <initializer_list>
#include <string>
#include <utility>
#include <vector>

namespace axom
{
namespace utilities
{
/*!
 * \brief Path class for performing basic path operations with user-selectable
 * delimiter characters
 */
class Path
{
public:
  /*!
     * \brief Constructs a path from a string
     *
     * \param [in] path The string representing the path
     * \param [in] delim The character that delimits the path
     * 
     * FIXME: Should this be marked \p explicit ?
     */
  Path(const std::string& path, const char delim = '/');

  /*!
     * \brief Default constructor for an empty path
     */
  Path() = default;

  /*!
     * \brief Constructs a path by concatenating a set of paths
     *
     * \param [in] paths The paths to concatenate/join
     * \param [in] delim The character to join the paths with
     * TODO: Should this be something easier-to-use, like the delimiter used
     * by the first Path object in the list?
     * 
     * \return The joined paths
     */
  static Path join(std::initializer_list<Path> paths, const char delim = '/');

  /*!
     * \brief Converts a path to a string
     * 
     * FIXME: Should this be marked \p explicit ?
     */
  operator std::string() const;

  /*!
     * \brief Returns a path representing the parent of the calling path
     */
  Path parent() const;

  /*!
     * \brief Returns the basename of the path (the last component)
     * \see basename(3)
     */
  std::string baseName() const;

  /*!
     * \brief Returns dirname of the path (all but the last component)
     * \see dirname(3)
     */
  std::string dirName() const;

  /*!
     * \brief Splits the path into the "parent" (dirname) and the last component (basename)
     * 
     * \return A pair {dirname, basename}
     */
  std::pair<std::string, std::string> split() const;

private:
  // The components (tokens) that make up the path
  std::vector<std::string> m_components;

  // The delimiter character
  char m_delim = '/';
};

}  // end namespace utilities
}  // end namespace axom

#endif  // CORE_PATH_UTILITIES_H_
