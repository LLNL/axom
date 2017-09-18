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
#include "conduit_relay.hpp"

#include "mpi.h"

using axom::spio::IOManager;
using axom::sidre::Group;
using axom::sidre::DataStore;
using axom::sidre::DataType;
using axom::sidre::View;


//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  SCR_Init();

  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  int num_ranks;
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  int num_output = num_ranks;

  /*
   * Create a DataStore and give it a small hierarchy of groups and views.
   *
   * The views are filled with repeatable nonsense data that will vary based
   * on rank.
   */
  DataStore * ds = new DataStore();

  Group * root = ds->getRoot();

  Group * flds = root->createGroup("fields");
  Group * flds2 = root->createGroup("fields2");

  Group * ga = flds->createGroup("a");
  Group * gb = flds2->createGroup("b");
  ga->createViewScalar<int>("i0", 101*my_rank);
  gb->createView("i1")->allocate(DataType::c_int(10));
  int * i1_vals = gb->getView("i1")->getData();

  for(int i=0 ; i<10 ; i++)
  {
    i1_vals[i] = (i+10) * (404-my_rank-i);
  }

  /*
   * Contents of the DataStore written to files with IOManager.
   */
  int num_files = num_output;
  IOManager writer(MPI_COMM_WORLD, true);

  writer.write(root, num_files, "out_spio_parallel_write_read", "sidre_hdf5");

  std::string root_name = "out_spio_parallel_write_read.root";

  MPI_Barrier(MPI_COMM_WORLD);

  int return_val = 0;

  /*
   * Create another DataStore that holds nothing but the root group.
   */
  DataStore * ds2 = new DataStore();

  /*
   * Read from the files that were written above.
   */
  IOManager reader(MPI_COMM_WORLD, true);

  reader.read(ds2->getRoot(), root_name, false, true);

  /*
   * Verify that the contents of ds2 match those written from ds.
   */
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

  View * view_i1_orig =
    ds->getRoot()->getGroup("fields2")->getGroup("b")->getView("i1");
  View * view_i1_restored =
    ds2->getRoot()->getGroup("fields2")->getGroup("b")->getView("i1");

  int num_elems = view_i1_orig->getNumElements();
  if (view_i1_restored->getNumElements() != num_elems) {
    return_val = 1;
  }
  else
  {
    int * i1_orig = view_i1_orig->getData();
    int * i1_restored = view_i1_restored->getData();

    for (int i = 0; i < num_elems; ++i) {
      if (return_val != 1) {
	if (i1_orig[i] != i1_restored[i]) {
	  return_val = 1;
	}
      }
    } 
  }

  delete ds;
  delete ds2;

  SCR_Finalize();
  MPI_Finalize();

  return return_val;
}

