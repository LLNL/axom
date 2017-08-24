/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */


/**************************************************************************
 *************************************************************************/

#include "mpi.h"

#include "slic/slic.hpp"
#include "slic/UnitTestLogger.hpp"

#include "sidre/Group.hpp"
#include "sidre/DataStore.hpp"
#include "spio/IOManager.hpp"

using axom::sidre::Group;
using axom::sidre::DataStore;
using axom::spio::IOManager;

/**************************************************************************
 * Subroutine:  main
 * Purpose   :
 *************************************************************************/

int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);
  SCR_Init();
  axom::slic::UnitTestLogger logger;

  SLIC_ERROR_IF(argc != 2,
      "Missing required command line argument. \n\t"
      << "Usage: spio_IORead <sidre_root_file>");

  DataStore * ds = new DataStore();
  SLIC_ASSERT(ds);
  Group * root = ds->getRoot();

  std::string root_file;
  if (argc == 2) {
    root_file = argv[1];
  } else {
    return 0;
  }

  IOManager reader(MPI_COMM_WORLD, true);
  reader.readWithSCR(root, root_file);

  delete ds;

  SCR_Finalize();
  MPI_Finalize();


  return 0;
}
