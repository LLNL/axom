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
#include "ioparallel/IOParallel.hpp"

using asctoolkit::sidre::DataGroup;
using asctoolkit::sidre::DataStore;
using asctoolkit::ioparallel::IOParallel;

/**************************************************************************
 * Subroutine:  main
 * Purpose   :
 *************************************************************************/

int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);

  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();
  DataGroup * flds = root->createGroup("fields");
  DataGroup * flds2 = root->createGroup("fields2");

  size_t num_files = 0;
  std::string file_base;
  if (argc == 3) {
    num_files = static_cast<size_t>(atoi(argv[1]));
    file_base = argv[2];
  } else {
    return 0;
  }

  std::vector<DataGroup *> groups;
  groups.push_back(root);

  IOParallel reader(MPI_COMM_WORLD, groups, num_files);
  reader.read(file_base, 0, "conduit");

  delete ds;

  MPI_Finalize();


  return 0;
}
