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

  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  int num_ranks;
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  int num_output = num_ranks / 2; 
  if (num_output == 0) {
    num_output = 1;
  }

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
  IOManager writer(MPI_COMM_WORLD);

  writer.write(root, num_files, "out_spio_preserve_write_read", "sidre_hdf5");

  std::string root_name = "out_spio_preserve_write_read.root";

  /*
   * Extra stuff to exercise preserve_contents option
   */
  MPI_Barrier(MPI_COMM_WORLD);
  DataStore * dsextra = new DataStore();
  Group * extra = dsextra->getRoot()->createGroup("extra");
  extra->createViewScalar<double>("dval", 1.1);
  Group * child = extra->createGroup("child");
  child->createViewScalar<int>("ival", 7);
  child->createViewString("word0", "hello");
  child->createViewString("word1", "world");

  writer.write(extra, num_files, "out_spio_extra", "sidre_hdf5");
  std::string extra_root = "out_spio_extra.root";
 
  MPI_Barrier(MPI_COMM_WORLD);

  int return_val = 0;

  /*
   * Create another DataStore that holds nothing but the root group.
   */
  DataStore * ds2 = new DataStore();

  /*
   * Read from the files that were written above.
   */
  IOManager reader(MPI_COMM_WORLD);

  reader.read(ds2->getRoot(), root_name);

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

  /*
   * Read in extra file while preserving contents
   */
  Group * extra_fields = ds2->getRoot()->getGroup("fields");
  reader.read(extra_fields, extra_root, true);

  /*
   * Test one of the pre-existing Views to show it's unchanged.
   */
  int testvalue_extra = extra_fields->getGroup("a")->getView("i0")->getData();
  if (testvalue_extra != testvalue2) {
    return_val = 1;
  }

  /*
   * Test the data from the extra file.
   */
  double dval = extra_fields->getView("dval")->getData();

  if (dval < 1.0000009 || dval > 1.1000001) {
    return_val = 1;
  }

  int ival = extra_fields->getView("child/ival")->getData();
  if (ival != 7) {
    return_val = 1;
  }

  if (std::string(extra_fields->getView("child/word0")->getString()) != "hello") {
    return_val = 1;
  }

  if (std::string(extra_fields->getView("child/word1")->getString()) != "world") {
    return_val = 1;
  }

  delete ds;
  delete ds2;
  delete dsextra;

  MPI_Finalize();

  return return_val;
}

