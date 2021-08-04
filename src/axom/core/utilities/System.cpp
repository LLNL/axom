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
#endif

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
  bufCharCount = INFO_BUFFER_SIZE;
  if(!GetComputerName(infoBuf, &bufCharCount))
  {
    hostName = std::string(infoBuf);
  }
#else
  char infoBuf[HOST_NAME_MAX];
  if(gethostname(infoBuf, HOST_NAME_MAX) == 0)
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
  bufCharCount = INFO_BUFFER_SIZE;
  if(GetUserName(infoBuf, &bufCharCount))
  {
    userName = std::string(infoBuf);
  }
#else
  char infoBuf[LOGIN_NAME_MAX];
  if(getlogin_r(infoBuf, LOGIN_NAME_MAX) == 0)
  {
    SLIC_INFO("~~~ getlogin_r succeeded!")
    userName = std::string(infoBuf);
  }
  else
  {
    SLIC_INFO("~~~ getlogin_r failed!")
  }
#endif
  return userName;
}

}  // end namespace utilities
}  // end namespace axom
