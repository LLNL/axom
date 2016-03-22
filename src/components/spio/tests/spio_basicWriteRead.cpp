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

  DataStore * ds = new DataStore();

  DataGroup * root = ds->getRoot();

  DataGroup * flds = root->createGroup("fields");
  DataGroup * flds2 = root->createGroup("fields2");

  DataGroup * ga = flds->createGroup("a");
  DataGroup * gb = flds2->createGroup("b");
  ga->createView("i0")->allocate(DataType::c_int());
  ga->getView("i0")->setScalar(101);
  gb->createView("i1")->allocate(DataType::c_int());
  gb->getView("i1")->setScalar(404);

  std::vector<DataGroup *> groups;
  groups.push_back(root);

  int num_files = 1;
  IOManager writer(MPI_COMM_WORLD, &(groups[0]), groups.size(), num_files);

  writer.write("out_spio_basic_write_read", 0, "conduit_hdf5");

  DataStore * ds2 = new DataStore();
  std::vector<DataGroup *> groups2;
  groups2.push_back(ds2->getRoot());

  IOManager writer2(MPI_COMM_WORLD, &(groups2[0]), groups2.size(), num_files);

  writer2.read("out_spio_basic_write_read", 0, "conduit_hdf5");

  int return_val = 0;
  if (!ds2->getRoot()->isEquivalentTo(root)) {
    return_val = 1; 
  }

  int testvalue =
    ds->getRoot()->getGroup("fields")->getGroup("a")->getView("i0")->getData();
  int testvalue2 =
    ds2->getRoot()->getGroup("fields")->getGroup("a")->getView("i0")->getData();

  if (testvalue != testvalue2) {
    return_val = 1;
  }

  testvalue =
    ds->getRoot()->getGroup("fields2")->getGroup("b")->getView("i1")->getData();
  testvalue2 =
    ds2->getRoot()->getGroup("fields2")->getGroup("b")->getView("i1")->getData();

  if (testvalue != testvalue2) {
    return_val = 1;
  }

  delete ds;
  delete ds2;

  MPI_Finalize();

  return return_val;
}

