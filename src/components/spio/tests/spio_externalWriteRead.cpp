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

  const int nvals = 10;
  int orig_vals1[nvals], orig_vals2[nvals];
  for(int i=0 ; i<10 ; i++)
  {
    orig_vals1[i] = (i+10) * (404-my_rank-i);
    orig_vals2[i] = (i+10) * (404-my_rank-i) + 20;
  }

  /*
   * Create a DataStore and give it a small hierarchy of groups and views.
   *
   * The views are filled with repeatable nonsense data that will vary based
   * on rank.
   */
  DataStore * ds1 = new DataStore();

  DataGroup * root1 = ds1->getRoot();

  DataGroup * flds = root1->createGroup("fields");
  DataGroup * flds2 = root1->createGroup("fields2");

  DataGroup * ga = flds->createGroup("a");
  DataGroup * gb = flds2->createGroup("b");
  ga->createView("external_array", asctoolkit::sidre::INT_ID, nvals, orig_vals1);
  gb->createView("external_undescribed")->setExternalDataPtr(orig_vals2);

  /*
   * Contents of the DataStore written to files with IOManager.
   */
  int num_files = num_output;
  IOManager writer(MPI_COMM_WORLD);

  writer.write(root1, num_files, "out_spio_external_write_read", "sidre_hdf5");

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

  int restored_vals1[nvals], restored_vals2[nvals];
  for (int i = 0; i < nvals; ++i) {
    restored_vals1[i] = -1;
    restored_vals2[i] = -1;
  }

  DataView * view1 = root2->getView("fields/a/external_array");
  view1->setExternalDataPtr(restored_vals1);

  DataView * view2 = root2->getView("fields2/b/external_undescribed");
  view2->setExternalDataPtr(restored_vals2);

  reader.loadExternalData(root2, "out_spio_external_write_read.root");


  /*
   * Verify that the contents of ds2 match those written from ds.
   */
  int return_val = 0;
  if (!ds2->getRoot()->isEquivalentTo(root1)) {
    return_val = 1;
  }

  if (view1->getNumElements() != nvals) {
    return_val = 1;
  } else {
      for (int i = 0; i < nvals; ++i) {
	  if (orig_vals1[i] != restored_vals1[i]) {
	      return_val = 1;
	      break;
	  }
      } 
  }

  /*
   * external_undescribed was not written to disk (since it is undescribed)
   * make sure it was not read in.
   */
  for (int i = 0; i < nvals; ++i) {
      if (-1 != restored_vals2[i]) {
	  return_val = 1;
	  break;
      }
  } 

  delete ds1;
  delete ds2;

  MPI_Finalize();

  return return_val;
}

