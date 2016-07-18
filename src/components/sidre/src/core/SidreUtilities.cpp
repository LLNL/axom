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
 * \brief   Implementation file for SidreUtilities.
 *
 ******************************************************************************
 */

// Standard C++ headers
#include <string>
#include <vector>

// Associated header file
#include "SidreUtilities.hpp"

// Other toolkit component headers
#include "slic/slic.hpp"

namespace asctoolkit
{
namespace sidre
{
namespace detail
{

/*
 * Splits a string given the provided delimiter into a vector of strings.
 * If the position of the first delimiter is known, it can be provided in the 'pos' parameter.
 */
std::vector<std::string> split(const std::string& s, char c,
                               std::string::size_type pos, bool keep_empty)
{
  std::vector<std::string> v;
  std::string::size_type i = 0;

  if (pos == std::string::npos)
  {
    pos = s.find(c);
  }
  while (pos != std::string::npos)
  {
    if (keep_empty || pos-i > 0)
    {
      v.push_back(s.substr(i, pos-i));
    }
    i = ++pos;
    pos = s.find(c, pos);

    if (pos == std::string::npos && (keep_empty || s.length() - i > 0))
    {
      v.push_back(s.substr(i, s.length()));
    }
  }

  return v;
}

/*
 * Checks if a specified character is found in a string, but a valid input
 * string is defined as only having the character in the internals of the
 * string, exclusive of the first and last position.
 *
 * It will return the position of the first instance found, or 0 if none
 * found.
 *
 * This function is used by DataGroup when checking for a valid path,
 * and the path delimiter is not considered valid when in the first
 * or last spot of the string.  See ATK-316 for more info.
 */
std::string::size_type find_exclusive( const std::string& s, char c)
{
  std::string::size_type pos = s.find(c);

  SLIC_ASSERT( pos != 0);
  SLIC_ASSERT( pos != s.size()-1 );

  if (pos == 0 || pos == s.size()-1)
  {
    pos = std::string::npos;
  }

  return pos;
}


} //end namespace detail

} //end namespace sidre

} // end namespace asctoolkit
