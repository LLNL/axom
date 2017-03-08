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

#include "common/FileUtilities.hpp"
#include "sidre/DataGroup.hpp"
#include "sidre/DataStore.hpp"
#include "spio/IOManager.hpp"

using axom::sidre::DataGroup;
using axom::sidre::DataStore;
using axom::sidre::DataType;
using axom::spio::IOManager;
using namespace axom::utilities;

/**************************************************************************
 * Subroutine:  main
 * Purpose   :
 *************************************************************************/

int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);
  axom::slic::UnitTestLogger logger;

  SLIC_ERROR_IF(argc != 3,
      "Missing command line arguments. \n\t"
      << "Usage: spio_IOWrite <num_files> <base_file_name>");

  size_t num_files = 0;
  std::string file_base;
  if (argc == 3) {
    num_files = static_cast<size_t>(atoi(argv[1]));
    file_base = argv[2];
  } else {
    return 0;
  }

  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  DataStore * ds = new DataStore();
  SLIC_ASSERT(ds);
  DataGroup * root = ds->getRoot();

  DataGroup * flds = root->createGroup("fields");
  DataGroup * flds2 = root->createGroup("fields2");

  DataGroup * ga = flds->createGroup("a");
  DataGroup * gb = flds2->createGroup("b");
  ga->createViewScalar<int>("i0", my_rank + 101);
  gb->createViewScalar<int>("i1", 4*my_rank*my_rank + 404);

  if (my_rank == 0) {
    std::string dir;
    filesystem::getDirName(dir, file_base);
    filesystem::makeDirsForPath(dir);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  IOManager writer(MPI_COMM_WORLD);
  writer.write(root, num_files, file_base, "sidre_hdf5");

  MPI_Barrier(MPI_COMM_WORLD);
  if (my_rank == 0) {
    DataGroup * extra = root->createGroup("extra");
    extra->createViewScalar<double>("dval", 1.1);
    DataGroup * child = extra->createGroup("child");
    child->createViewScalar<int>("ival", 7);
    child->createViewString("word0", "hello");
    child->createViewString("word1", "world");

    std::string root_name = file_base + ".root";
    writer.writeGroupToRootFile(extra, root_name);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  delete ds;

  MPI_Finalize();


  return 0;
}
