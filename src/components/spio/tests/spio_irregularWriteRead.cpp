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

#include "mpi.h"

using axom::spio::IOManager;
using axom::sidre::DataGroup;
using axom::sidre::DataStore;
using axom::sidre::DataType;
using axom::sidre::DataView;


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
  DataStore * ds1 = new DataStore();

  DataGroup * root1 = ds1->getRoot();

  int num_fields = my_rank + 2;

  for (int f = 0; f < num_fields; ++f) {
    std::ostringstream ostream;
    ostream << "fields" << f;
    DataGroup * flds = root1->createGroup(ostream.str());

    int num_subgroups = ((f+my_rank)%3) + 1;
    for (int g = 0; g < num_subgroups; ++g) {
      std::ostringstream gstream;
      gstream << "subgroup" << g;
      DataGroup * sg = flds->createGroup(gstream.str());

      std::ostringstream vstream;
      vstream << "view" << g;
      if (g % 2) {
        sg->createView(vstream.str())->allocate(DataType::c_int(10+my_rank));
        int * vals = sg->getView(vstream.str())->getData();

        for(int i=0 ; i<10+my_rank ; i++)
        {
          vals[i] = (i+10) * (404-my_rank-i-g-f);
        }

      } else {
        sg->createViewScalar<int>(vstream.str(), 101*my_rank*(f+g+1));
      } 
    } 
  } 


  /*
   * Contents of the DataStore written to files with IOManager.
   */
  int num_files = num_output;
  IOManager writer(MPI_COMM_WORLD);

  writer.write(root1, num_files, "out_spio_irregular_write_read", "sidre_hdf5");

  /*
   * Create another DataStore that holds nothing but the root group.
   */
  DataStore * ds2 = new DataStore();

  /*
   * Read from the files that were written above.
   */
  IOManager reader(MPI_COMM_WORLD);

  reader.read(ds2->getRoot(), "out_spio_irregular_write_read.root");


  /*
   * Verify that the contents of ds2 match those written from ds.
   */
  int return_val = 0;
  if (!ds2->getRoot()->isEquivalentTo(root1)) {
    return_val = 1; 
  }

  for (int f = 0; f < num_fields; ++f) {
    std::ostringstream ostream;
    ostream << "fields" << f;
    DataGroup * flds1 = ds1->getRoot()->getGroup(ostream.str());
    DataGroup * flds2 = ds2->getRoot()->getGroup(ostream.str());

    int num_subgroups = ((f+my_rank)%3) + 1;
    for (int g = 0; g < num_subgroups; ++g) {
      std::ostringstream gstream;
      gstream << "subgroup" << g;
      DataGroup * sg1 = flds1->getGroup(gstream.str());
      DataGroup * sg2 = flds2->getGroup(gstream.str());

      std::ostringstream vstream;
      vstream << "view" << g;
      if (g % 2) {

        DataView * view_orig = sg1->getView(vstream.str());
        DataView * view_restored = sg2->getView(vstream.str());

        int num_elems = view_orig->getNumElements();
        if (view_restored->getNumElements() != num_elems) {
          return_val = 1; 
        }
        else
        {
	  int * vals_orig = view_orig->getData();
	  int * vals_restored = view_restored->getData();

	  for (int i = 0; i < num_elems; ++i) {
	    if (return_val != 1) {
	      if (vals_orig[i] != vals_restored[i]) {
		return_val = 1;
	      }
	    }
	  }
	}
      } else {
        int testvalue1 = sg1->getView(vstream.str())->getData();
        int testvalue2 = sg2->getView(vstream.str())->getData();

        if (testvalue1 != testvalue2) {
          return_val = 1;
        }
      }
    }
  } 

  delete ds1;
  delete ds2;

  MPI_Finalize();

  return return_val;
}

