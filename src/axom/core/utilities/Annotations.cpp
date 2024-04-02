#include "Annotations.hpp"

#ifdef AXOM_USE_MPI
  #include <mpi.h>
#endif

namespace axom
{
namespace utilities
{
namespace annotations
{
static bool adiak_initialized = false;

void initialize()
{
#ifdef AXOM_USE_ADIAK
  if(adiak_initialized)
  {
    return;
  }

  #ifdef AXOM_USE_MPI
  MPI_Comm communicator = MPI_COMM_WORLD;
  adiak::init((void*)&communicator);
  #else
  adiak::init(nullptr);
  #endif

  adiak::launchdate();
  adiak::executable();
  adiak::cmdline();
  adiak::clustername();
  adiak::jobsize();
  adiak::walltime();
  adiak::cputime();
  adiak::systime();

  adiak_initialized = true;
#endif
}

void finalize()
{
  if(adiak_initialized)
  {
    adiak::fini();
  }
  adiak_initialized = false;
}

}  // namespace annotations
}  // namespace utilities
}  // namespace axom