// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file SidreUtilities.cpp
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

// Other axom headers
#include "axom/slic/interface/slic.hpp"

/// @cond INCLUDE_DETAIL
namespace axom
{
namespace sidre
{
namespace detail
{

/*
 *************************************************************************
 *
 * Split a string s (such as a pathname) into substrings delimited by c
 *
 *************************************************************************
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
 *************************************************************************
 *
 * Report the position of the first occurence of c in s.  First or last
 * or not found is reported with std::string::npos.  In addition, c being
 * first or last is considered an error.
 *
 *************************************************************************
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

} // end namespace axom
/// @endcond
