// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

//-----------------------------------------------------------------------------
//
// file: utils_locale.hpp
// Checks locales and prints out related info
//
//-----------------------------------------------------------------------------

#include "axom/core/utilities/StringUtilities.hpp"

#include "axom/fmt.hpp"

#include <locale>
#include <string>
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <array>
#include <set>
#include <algorithm>

#include "gtest/gtest.h"

namespace
{

#ifdef __linux__
/**
 * \brief Utility function to capture output of system call
 * \note Adapted from https://stackoverflow.com/a/478960
*/
std::string exec(const char* cmd)
{
  std::array<char, 128> buffer;
  std::string result;

  std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
  if(!pipe)
  {
    throw std::runtime_error(
      axom::fmt::format("popen() failed for command '{}'", cmd));
  }

  while(fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr)
  {
    result += buffer.data();
  }
  return result;
}
#endif  // __linux__

}  // end anonymous namespace

//-----------------------------------------------------------------------------
TEST(utils_locale, default_locale)
{
  std::locale loc(std::locale(), new std::ctype<char>);
  std::cout << "The default locale is " << std::locale().name() << '\n'
            << "The user's locale is " << std::locale("").name() << '\n'
            << "Nameless locale is " << loc.name() << '\n';
}

//-----------------------------------------------------------------------------
TEST(utils_locale, has_en_US_locale)
{
  try
  {
    std::locale loc("en_US.UTF-8");
    std::cout << "Name after setting locale to 'en_US.UTF-8': " << loc.name()
              << std::endl;
  }
  catch(std::runtime_error&)
  {
    std::cout << "Could not set locale to 'en_US.UTF-8'" << std::endl;
  }
}

TEST(utils_locale, try_catch)
{
  std::locale loc;  // initialized to locale::classic()

  try
  {
    loc = std::locale("en_US.UTF8");
  }
  catch(std::runtime_error&)
  {
    loc = std::locale(loc, "", std::locale::ctype);
  }

  std::cout << "The selected locale after attempting to use 'en_US.UTF-8' "
               "with try/catch is: "
            << loc.name() << std::endl;
}

#ifdef __linux__
TEST(utils_locale, enumerate_locales_linux)
{
  std::vector<std::string> locale_list;
  std::vector<std::string> available_locales;

  // get list of locales by running "locale -a" system call
  try
  {
    std::string output = ::exec("locale -a");
    locale_list = axom::utilities::string::split(output, '\n');
  }
  catch(const std::runtime_error& e)
  {
    std::cerr << e.what() << '\n';
  }

  // sort and unique-ify the list
  std::sort(locale_list.begin(), locale_list.end());
  locale_list.erase(std::unique(locale_list.begin(), locale_list.end()),
                    locale_list.end());

  // check if each item can successfully be used to create a locale
  for(const auto& loc_str : locale_list)
  {
    try
    {
      auto loc = std::locale(loc_str);
      available_locales.push_back(loc.name());
    }
    catch(std::runtime_error&)
    {
      // could not instantiate this locale
    }
  }

  std::cout << axom::fmt::format("Available locales on this system: {}",
                                 axom::fmt::join(available_locales, ", "))
            << std::endl;
}
#endif  // __linux__