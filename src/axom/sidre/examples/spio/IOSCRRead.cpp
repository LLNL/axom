// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**************************************************************************
 *************************************************************************/

#include "mpi.h"

#include "axom/slic.hpp"
#include "axom/sidre.hpp"

#ifdef AXOM_USE_SCR
#include "scr.h"
#endif

using axom::sidre::Group;
using axom::sidre::DataStore;
using axom::sidre::IOManager;

/**************************************************************************
 * Subroutine:  main
 * Purpose   :
 *************************************************************************/

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  SCR_Init();
  axom::slic::UnitTestLogger logger;

  SLIC_ERROR_IF(argc != 2,
                "Missing required command line argument. \n\t"
                << "Usage: spio_IORead <sidre_root_file>");

  DataStore* ds = new DataStore();
  SLIC_ASSERT(ds);
  Group* root = ds->getRoot();

  std::string root_file;
  if (argc == 2)
  {
    root_file = argv[1];
  }
  else
  {
    return 0;
  }

  IOManager reader(MPI_COMM_WORLD, true);
  reader.read(root, root_file, false, true);

  delete ds;

  SCR_Finalize();
  MPI_Finalize();


  return 0;
}
