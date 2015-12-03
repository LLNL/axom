#include "common/Utilities.hpp"

#include <cstdlib>

#ifdef USE_MPI
#include <mpi.h>
#endif

namespace asctoolkit
{
namespace utilities
{

void processAbort()
{
#ifndef USE_MPI
  exit( -1 );
#else
  int mpi = 0;
  MPI_Initialized( &mpi );
  if ( mpi ) {

    MPI_Abort( MPI_COMM_WORLD, -1 );

  } else {

    exit( -1 );

  }
#endif
}

}
}
