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

  writer.write(root, num_files, "out_spio_parallel_write_read", "sidre_hdf5");

  std::string root_name = "out_spio_parallel_write_read.root";

  /*
   * Extra stuff to exercise writeGroupToRootFile
   */
  MPI_Barrier(MPI_COMM_WORLD);
  if (my_rank == 0) { 
    DataStore * dsextra = new DataStore();
    DataGroup * extra = dsextra->getRoot()->createGroup("extra");
    extra->createViewScalar<double>("dval", 1.1);
    DataGroup * child = extra->createGroup("child");
    child->createViewScalar<int>("ival", 7);
    child->createViewString("word0", "hello");
    child->createViewString("word1", "world");

    writer.writeGroupToRootFile(extra, root_name);

    DataGroup * path_test = dsextra->getRoot()->createGroup("path_test");

    path_test->createViewScalar<int>("path_val", 9);
    path_test->createViewString("word2", "again");

    writer.writeGroupToRootFileAtPath(path_test, root_name, "extra/child");

    DataView * view_test = dsextra->getRoot()->createViewString("word3", "new_view");

    writer.writeViewToRootFileAtPath(view_test,
                                     root_name,
                                     "extra/child/path_test");

    delete dsextra;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  int return_val = 0;
  /*
   * Read the root file on rank 1, unless this is a serial run.
   */
  if (my_rank == 1 || num_ranks == 1) {

    conduit::Node n;
    conduit::relay::io::load("hdf5", root_name + ":extra", n);

    double dval = n["dval"].to_double();

    if (dval < 1.0000009 || dval > 1.1000001) {
      return_val = 1;
    }

    if (n["child"]["ival"].to_int() != 7) {
      return_val = 1;
    }

    if (n["child"]["word0"].as_string() != "hello") {
      return_val = 1;
    }

    if (n["child"]["word1"].as_string() != "world") {
      return_val = 1;
    }

    if (n["child"]["path_test"]["path_val"].to_int() != 9) {
      return_val = 1;
    }

    if (n["child"]["path_test"]["word2"].as_string() != "again") {
      return_val = 1;
    }

    if (n["child"]["path_test"]["word3"].as_string() != "new_view") {
      return_val = 1;
    }

  }

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

  DataView * view_i1_orig =
    ds->getRoot()->getGroup("fields2")->getGroup("b")->getView("i1");
  DataView * view_i1_restored =
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

  MPI_Finalize();

  return return_val;
}

