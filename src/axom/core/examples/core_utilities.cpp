// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! \file core_utilities.cpp
 *  \brief This example code is a demonstration of the Axom Core utilites.
 *
 *  This file shows how to use Core::utilities
 */

 /* This example code contains snippets used in the Primal Sphinx documentation.
  * They begin and end with comments such as
  *
  * timer_start
  * timer_end
  *
  * each prepended with an underscore.
  */


#ifdef WIN32
#include "windows.h"
void sleep(int numSeconds)
{
  int numMilliSecs = numSeconds * 1000;
  Sleep(numMilliSecs);
}
#else
#include <unistd.h> // for sleep()
#endif

// C/C++ includes
#include <iostream>
#include <vector>

// Axom includes
// #include "axom/core.hpp"
#include "axom/core/utilities/AnnotationMacros.hpp"
#include "axom/core/utilities/FileUtilities.hpp"
#include "axom/core/utilities/StringUtilities.hpp"
#include "axom/core/utilities/Timer.hpp"
// #include "axom/core/utilities/Utilities.hpp"
#include "axom/core/utilities/About.hpp"

// _using_start
// _using_end

// _fs_string_start
void demoFileSystemAndString(const char* argv0)
{
  using namespace axom::utilities;

  // Get the current directory
  std::string cwd = filesystem::getCWD();

  // Split it on file separator.
#if WIN32
  const char pathsep = '\\';
#else
  const char pathsep = '/';
#endif
  std::vector<std::string> cmp;
  string::split(cmp, cwd, pathsep);

  // Count how many start with "ax" or end with "exe"
  // (we could also use std::count_if)
  int matchcount = 0;
  std::string prefix{ "ax" };
  std::string suffix{ "exe" };
  for (int i = 0; i < cmp.size(); ++i)
  {
    if (string::startsWith(cmp[i], prefix) || 
        string::endsWith(cmp[i], suffix))
    {
      matchcount += 1;
    }
  }
  std::cout << "Found " << matchcount << " path components starting with " << 
    prefix << " or ending with " << suffix << "." << std::endl;

  // Append "hello.txt"
  std::string hellofile = 
    filesystem::joinPath(cwd, "hello.txt", std::string{pathsep});

  // Does this exist?
  std::cout << "The file \"hello.txt\" ";
  if (filesystem::pathExists(hellofile))
  {
    std::cout << "exists ";
  }
  else
  {
    std::cout << "DOES NOT exist ";
  }
  std::cout << "in the current working directory." << std::endl;

  // Does argv0 exist?
  if (filesystem::pathExists(argv0))
  {
    std::cout << argv0 << " exists ";
  }
  else
  {
    std::cout << argv0 << " DOES NOT exist ";
  }
  std::cout << "in the current working directory." << std::endl;

  // sleep for a second
  sleep(1);
}
// _fs_string_end

int main(int argc, char** argv)
{
  // _timer_start
  axom::utilities::Timer t;

  t.start();

  demoFileSystemAndString(argv[0]);

  t.stop();

  std::cout << "The tests took " << t.elapsedTimeInMilliSec() <<
    " ms." << std::endl;
  // _timer_end
}
