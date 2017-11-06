/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include "gtest/gtest.h"

#include "sidre/sidre.hpp"
#include "spio/IOManager.hpp"

using axom::sidre::Group;
using axom::sidre::DataStore;
using axom::spio::IOManager;

//------------------------------------------------------------------------------

TEST(spio_serial, basic_writeread)
{
  DataStore * ds1 = new DataStore();

  Group * root1 = ds1->getRoot();

  Group * flds1 = root1->createGroup("fields");
  Group * flds2 = root1->createGroup("fields2");

  Group * ga = flds1->createGroup("a");
  Group * gb = flds2->createGroup("b");
  ga->createViewScalar<int>("i0", 101);
  gb->createViewScalar<int>("i1", 404);

  int num_files = 1;
  IOManager writer(MPI_COMM_WORLD);

  writer.write(root1, num_files, "out_spio_basic_write_read", "sidre_hdf5");

  DataStore * ds2 = new DataStore();

  IOManager reader(MPI_COMM_WORLD);

  reader.read(ds2->getRoot(), "out_spio_basic_write_read.root");

  EXPECT_TRUE(ds2->getRoot()->isEquivalentTo(root1));

  int testvalue1 =
    ds1->getRoot()->getGroup("fields")->getGroup("a")->getView("i0")->getData();
  int testvalue2 =
    ds2->getRoot()->getGroup("fields")->getGroup("a")->getView("i0")->getData();

  EXPECT_EQ(testvalue1,testvalue2);

  testvalue1 =
    ds1->getRoot()->getGroup("fields2")->getGroup("b")->getView("i1")->getData();
  testvalue2 =
    ds2->getRoot()->getGroup("fields2")->getGroup("b")->getView("i1")->getData();

  EXPECT_EQ(testvalue1,testvalue2);

  delete ds1;
  delete ds2;
}

const std::string PROTOCOL = "sidre_hdf5";
const std::string ROOT_EXT = ".root";

//------------------------------------------------------------------------------
TEST(spio_serial, write_read_write)
{
  int my_rank, num_ranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  int num_files = std::max( num_ranks / 2, 1);
  std::stringstream sstr;
  sstr << "out_spio_WRW_" << num_ranks;
  std::string filename = sstr.str();

  // Initialize a datastore and dump to disk
  DataStore * ds = new DataStore();
  ds->getRoot()->createViewScalar("grp/i",2);
  ds->getRoot()->createViewScalar("grp/f",3.0);
  IOManager writer_a(MPI_COMM_WORLD);
  writer_a.write(ds->getRoot(), num_files, filename, PROTOCOL); 

  // Create another DataStore to read into.
  DataStore ds_r;
  IOManager reader(MPI_COMM_WORLD);
  reader.read(ds_r.getRoot(), filename + ROOT_EXT);

  // Dump this datastore to disk.
  // Regression: This used to produce the following HDF5 error:
  //  HDF5-DIAG: Error detected in HDF5 (1.8.16) thread 0:
  //    #000: H5F.c line 522 in H5Fcreate(): unable to create file
  //      major: File accessibility
  //      minor: Unable to open file
  //    #001: H5Fint.c line 1024 in H5F_open(): unable to truncate a file which is already open
  //      major: File accessibility
  //      minor: Unable to open file
  IOManager writer_b(MPI_COMM_WORLD);
  writer_b.write(ds_r.getRoot(), num_files, filename, PROTOCOL);

}



//----------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char * argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  MPI_Init(&argc, &argv);
  result = RUN_ALL_TESTS();
  MPI_Finalize();

  return result;
}
