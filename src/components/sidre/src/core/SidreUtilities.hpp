/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

#ifndef SIDREUTILITIES_HPP_
#define SIDREUTILITIES_HPP_

// Standard C++ headers
#include <string>
#include <vector>

namespace asctoolkit
{
namespace sidre
{
namespace detail
{
  std::vector<std::string> split(const std::string& s, char c, 
                                 std::string::size_type pos);

  std::string::size_type find_exclusive( const std::string& s, char c);

} // end namespace detail

} // end namespace sidre

} // end namespace asctoolkit

#endif /* UTILITIES_HPP_ */
