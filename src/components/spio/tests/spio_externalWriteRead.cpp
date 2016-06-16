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

  int nvals = 10;
  int orig_vals[nvals];
  for(int i=0 ; i<10 ; i++)
  {
    orig_vals[i] = (i+10) * (404-my_rank-i);
  }

  /*
   * Create a DataStore and give it a small hierarchy of groups and views.
   *
   * The views are filled with repeatable nonsense data that will vary based
   * on rank.
   */
  DataStore * ds = new DataStore();

  DataGroup * root = ds->getRoot();

  DataGroup * flds = root->createGroup("fields");
  DataGroup * flds2 = root->createGroup("fields2");

  DataGroup * ga = flds->createGroup("a");
  DataGroup * gb = flds2->createGroup("b");
  ga->createView("external_array", asctoolkit::sidre::INT_ID, nvals, orig_vals);
  gb->createView("external_undescribed")->setExternalDataPtr(orig_vals);

  /*
   * Contents of the DataStore written to files with IOManager.
   */
  int num_files = num_output;
  IOManager writer(MPI_COMM_WORLD);

  writer.write(root, num_files, "out_spio_external_write_read", "conduit_hdf5");

  /*
   * Create another DataStore than holds nothing but the root group.
   */
  DataStore * ds2 = new DataStore();
  DataGroup * root2 = ds2->getRoot();

  /*
   * Read from the files that were written above.
   */
  IOManager reader(MPI_COMM_WORLD);

  reader.read(root2, "out_spio_external_write_read.root");

  DataView * view1 = root2->getGroup("fields")->getGroup("a")->getView("external_array");

  int restored_vals[nvals];
  for (int i = 0; i < nvals; ++i) {
    restored_vals[i] = -1;
  }

  view1->setExternalDataPtr(restored_vals);

  DataView * view2 = root2->getGroup("fields2")->getGroup("b")->getView("external_undescribed");
  view2->setExternalDataPtr(restored_vals);

  reader.loadExternalData(root2, "out_spio_external_write_read.root"); 


  /*
   * Verify that the contents of ds2 match those written from ds.
   */
  int return_val = 0;
  if (!ds2->getRoot()->isEquivalentTo(root)) {
    return_val = 1;
  }

  if (view1->getNumElements() != nvals) {
    return_val = 1;
  }

  for (int i = 0; i < nvals; ++i) {
    if (return_val != 1) {
      if (orig_vals[i] != restored_vals[i]) {
        return_val = 1;
      }
    }
  } 

  delete ds;
  delete ds2;

  MPI_Finalize();

  return return_val;
}

