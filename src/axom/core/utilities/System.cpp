// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/utilities/System.hpp"

#ifdef WIN32
  #include <windows.h>
  #include <tchar.h>
#else
  #include <unistd.h>
  #include <limits.h>
  #include <pwd.h>
#endif

#include <iostream>

namespace axom
{
namespace utilities
{
#define INFO_BUFFER_SIZE 32767

std::string getHostName()
{
  std::string hostName = "";
#ifdef WIN32
  TCHAR infoBuf[INFO_BUFFER_SIZE];
  DWORD bufCharCount = INFO_BUFFER_SIZE;
  if(!GetComputerName(infoBuf, &bufCharCount))
  {
    hostName = std::string(infoBuf);
  }
#else
  char infoBuf[INFO_BUFFER_SIZE];
  if(gethostname(infoBuf, INFO_BUFFER_SIZE) == 0)
  {
    hostName = std::string(infoBuf);
  }
#endif
  return hostName;
}

std::string getUserName()
{
  std::string userName = "";
#ifdef WIN32
  TCHAR infoBuf[INFO_BUFFER_SIZE];
  DWORD bufCharCount = INFO_BUFFER_SIZE;
  if(GetUserName(infoBuf, &bufCharCount))
  {
    userName = std::string(infoBuf);
  }
#else
  char infoBuf[INFO_BUFFER_SIZE];
  if(getlogin_r(infoBuf, INFO_BUFFER_SIZE) == 0)
  {
    userName = std::string(infoBuf);
  }
  else
  {
    // fallback on getpwuid if getlogin_r fails
    struct passwd *pwd = getpwuid(getuid());
    if(pwd)
    {
      userName = std::string(pwd->pw_name);
    }
  }
#endif
  return userName;
}

}  // end namespace utilities
}  // end namespace axom
