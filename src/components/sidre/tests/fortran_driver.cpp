//----------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using asctoolkit::slic::UnitTestLogger;

extern "C" int fortran_test();

int main()
{
  int result = 0;

  UnitTestLogger logger;   // create & initialize test logger,
  // finalized when exiting main scope

  result = fortran_test();

  return result;
}
