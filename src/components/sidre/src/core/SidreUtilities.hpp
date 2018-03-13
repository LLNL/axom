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

#ifndef SIDREUTILITIES_HPP_
#define SIDREUTILITIES_HPP_

// Standard C++ headers
#include <string>
#include <vector>

namespace axom
{
namespace sidre
{
namespace detail
{

/*! \brief Splits a string using the given delimiter into a vector of strings.
 *
 *  \param s The string to split
 *  \param c Delimiter character
 *  \param pos Starting position (std::string::npos to start from beginning)
 *  \param keep_empty Include empty substrings in the result (default false)
 */
std::vector<std::string> split(const std::string& s, char c,
                               std::string::size_type pos,
                               bool keep_empty = false);

/*! \brief Report the first position of a character in a string,
 *  with first or last being an error.
 *
 *  \param s The string to scan
 *  \param c The character to scan for
 *
 * This function reports the position of the first instance of c found in s.
 * If c is not found, or if c is in the first or last position of s, the
 * function returns std::string::npos.
 *
 * This function is used by Group when checking for a valid path,
 * and the path delimiter is not considered valid when in the first
 * or last spot of the string.
 */
std::string::size_type find_exclusive( const std::string& s, char c);

} // end namespace detail

} // end namespace sidre

} // end namespace axom

#endif /* UTILITIES_HPP_ */
