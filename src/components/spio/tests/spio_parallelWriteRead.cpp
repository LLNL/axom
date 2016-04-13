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
using asctoolkit::sidre::DataView;


//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);

  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  int num_ranks;
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  int num_output = num_ranks / 2; 
  if (num_output == 0) {
    num_output = 1;
  }

  DataStore * ds = new DataStore();

  DataGroup * root = ds->getRoot();

  DataGroup * flds = root->createGroup("fields");
  DataGroup * flds2 = root->createGroup("fields2");

  DataGroup * ga = flds->createGroup("a");
  DataGroup * gb = flds2->createGroup("b");
  ga->createView("i0")->allocate(DataType::c_int());
  ga->getView("i0")->setScalar(101*my_rank);
  gb->createView("i1")->allocate(DataType::c_int(10));
  int * i1_vals = gb->getView("i1")->getData();

  for(int i=0 ; i<10 ; i++)
  {
    i1_vals[i] = (i+10) * (404-my_rank-i);
  }

  std::vector<DataGroup *> groups;
  groups.push_back(root);

  int num_files = num_output;
  IOManager writer(MPI_COMM_WORLD, &(groups[0]), groups.size(), num_files);

  writer.write("out_spio_parallel_write_read", 0, "conduit_hdf5");

  DataStore * ds2 = new DataStore();
  std::vector<DataGroup *> groups2;
  groups2.push_back(ds2->getRoot());

  IOManager writer2(MPI_COMM_WORLD, &(groups2[0]), groups2.size(), num_files);

  writer2.read("out_spio_parallel_write_read", 0, "conduit_hdf5");

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

  DataView * view_i1_orig =
    ds->getRoot()->getGroup("fields2")->getGroup("b")->getView("i1");
  DataView * view_i1_restored =
    ds2->getRoot()->getGroup("fields2")->getGroup("b")->getView("i1");

  int num_elems = view_i1_orig->getNumElements();
  if (view_i1_restored->getNumElements() != num_elems) {
    return_val = 1;
  }

  int * i1_orig = view_i1_orig->getData();
  int * i1_restored = view_i1_restored->getData();

  for (int i = 0; i < num_elems; ++i) {
    if (return_val != 1) {
      if (i1_orig[i] != i1_restored[i]) {
        return_val = 1;
      }
    }
  } 

  delete ds;
  delete ds2;

  MPI_Finalize();

  return return_val;
}

