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

#include "sidre/DataGroup.hpp"
#include "sidre/DataStore.hpp"
#include "spio/IOManager.hpp"

using asctoolkit::sidre::DataGroup;
using asctoolkit::sidre::DataStore;
using asctoolkit::spio::IOManager;

/**************************************************************************
 * Subroutine:  main
 * Purpose   :
 *************************************************************************/

int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);

  SLIC_ASSERT(argc == 2);

  DataStore * ds = new DataStore();
  SLIC_ASSERT(ds);
  DataGroup * root = ds->getRoot();

  std::string root_file;
  if (argc == 2) {
    root_file = argv[1];
  } else {
    return 0;
  }

  IOManager reader(MPI_COMM_WORLD);
  reader.read(root, root_file);

  delete ds;

  MPI_Finalize();


  return 0;
}
