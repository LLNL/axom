#ifndef MESHAPI_FILE_UTILITIES_H_
#define MESHAPI_FILE_UTILITIES_H_

#include <string>

#include <cstdio>  /* defines FILENAME_MAX */
#ifdef WINDOWS
    #include <direct.h>
    #define GetCurrentDir _getcwd
#else
    #include <unistd.h>
    #define GetCurrentDir getcwd
#endif


namespace asctoolkit {
namespace meshapi {
namespace util {

/**
 * \brief Helper function to print out the current working directory within the file system
 */
  std::string getCWD()
  {
    char cCurrentPath[FILENAME_MAX];

    if (!GetCurrentDir(cCurrentPath, FILENAME_MAX))
    {
      std::abort();
    }

    return std::string(cCurrentPath);
  }


} // end namespace util
} // end namespace meshapi
} // end namespace asctoolkit

#endif //  MESHAPI_FILE_UTILITIES_H_
