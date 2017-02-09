/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

#include "spio/IOManager.hpp"
#include "sidre/sidre.hpp"

using asctoolkit::spio::IOManager;
using asctoolkit::sidre::DataGroup;
using asctoolkit::sidre::DataStore;
using asctoolkit::sidre::DataType;


//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);

  DataStore * ds1 = new DataStore();

  DataGroup * root1 = ds1->getRoot();

  DataGroup * flds1 = root1->createGroup("fields");
  DataGroup * flds2 = root1->createGroup("fields2");

  DataGroup * ga = flds1->createGroup("a");
  DataGroup * gb = flds2->createGroup("b");
  ga->createViewScalar<int>("i0", 101);
  gb->createViewScalar<int>("i1", 404);

  int num_files = 1;
  IOManager writer(MPI_COMM_WORLD);

  writer.write(root1, num_files, "out_spio_basic_write_read", "sidre_hdf5");

  DataStore * ds2 = new DataStore();

  IOManager reader(MPI_COMM_WORLD);

  reader.read(ds2->getRoot(), "out_spio_basic_write_read.hdf5.root");
  
  int return_val = 0;
  if (!ds2->getRoot()->isEquivalentTo(root1)) {
    return_val = 1; 
  }

  int testvalue1 =
    ds1->getRoot()->getGroup("fields")->getGroup("a")->getView("i0")->getData();
  int testvalue2 =
    ds2->getRoot()->getGroup("fields")->getGroup("a")->getView("i0")->getData();

  if (testvalue1 != testvalue2) {
    return_val = 1;
  }

  testvalue1 =
    ds1->getRoot()->getGroup("fields2")->getGroup("b")->getView("i1")->getData();
  testvalue2 =
    ds2->getRoot()->getGroup("fields2")->getGroup("b")->getView("i1")->getData();

  if (testvalue1 != testvalue2) {
    return_val = 1;
  }

  delete ds1;
  delete ds2;

  MPI_Finalize();

  return return_val;
}

