// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"      // for AXOM_USE_HDF5
#include "axom/core/Types.hpp"  // for common::int64
#include "axom/sidre/core/sidre.hpp"
#include "axom/sidre/spio/IOManager.hpp"
#include <list>

#include "mpi.h"

using axom::sidre::DataStore;
using axom::sidre::Group;
using axom::sidre::IOManager;

using axom::int64;

//------------------------------------------------------------------------------

TEST(spio_serial, basic_writeread)
{
  DataStore* ds1 = new DataStore();

  Group* root1 = ds1->getRoot();

  Group* flds1 = root1->createGroup("fields");
  Group* flds2 = root1->createGroup("fields2");

  Group* ga = flds1->createGroup("a");
  Group* gb = flds2->createGroup("b");

  // Note: use 64-bit integers since that is the native type for json
  ga->createViewScalar<int64>("i0", 101);
  gb->createViewScalar<int64>("i1", 404);

  int num_files = 1;
  IOManager writer(MPI_COMM_WORLD);

  const std::string file_name = "out_spio_basic_write_read";

  writer.write(root1, num_files, file_name, PROTOCOL);

  DataStore* ds2 = new DataStore();

  IOManager reader(MPI_COMM_WORLD);

  reader.read(ds2->getRoot(), file_name + ROOT_EXT);

  EXPECT_TRUE(ds2->getRoot()->isEquivalentTo(root1));

  int64 testvalue1 =
    ds1->getRoot()->getGroup("fields")->getGroup("a")->getView("i0")->getData();
  int64 testvalue2 =
    ds2->getRoot()->getGroup("fields")->getGroup("a")->getView("i0")->getData();

  EXPECT_EQ(testvalue1, testvalue2);

  testvalue1 =
    ds1->getRoot()->getGroup("fields2")->getGroup("b")->getView("i1")->getData();
  testvalue2 =
    ds2->getRoot()->getGroup("fields2")->getGroup("b")->getView("i1")->getData();

  EXPECT_EQ(testvalue1, testvalue2);

  delete ds1;
  delete ds2;
}

//------------------------------------------------------------------------------

TEST(spio_serial, basic_writeread_protocols)
{
  std::list<std::string> protocols;
#ifdef AXOM_USE_HDF5
  protocols.push_back("conduit_hdf5");
#endif
  protocols.push_back("conduit_bin");
  protocols.push_back("conduit_json");
  protocols.push_back("json");
  protocols.push_back("sidre_conduit_json");
  protocols.push_back("sidre_json");

  for(std::list<std::string>::const_iterator itr = protocols.begin();
      itr != protocols.end();
      ++itr)
  {
    const std::string& protocol = *itr;

    DataStore* ds1 = new DataStore();

    Group* root1 = ds1->getRoot();

    Group* flds1 = root1->createGroup("fields");
    Group* flds2 = root1->createGroup("fields2");

    Group* ga = flds1->createGroup("a");
    Group* gb = flds2->createGroup("b");

    // Note: use 64-bit integers since that is the native type for json
    ga->createViewScalar<int64>("i0", 101);
    gb->createViewScalar<int64>("i1", 404);

    int num_files = 1;
    IOManager writer(MPI_COMM_WORLD);

    const std::string file_name = "out_spio_basic_write_read" + protocol;

    writer.write(root1, num_files, file_name, protocol);

    DataStore* ds2 = new DataStore();

    IOManager reader(MPI_COMM_WORLD);

    reader.read(ds2->getRoot(), file_name + ROOT_EXT);

    EXPECT_TRUE(ds2->getRoot()->isEquivalentTo(root1));

    int64 testvalue1 =
      ds1->getRoot()->getGroup("fields")->getGroup("a")->getView("i0")->getData();
    int64 testvalue2 =
      ds2->getRoot()->getGroup("fields")->getGroup("a")->getView("i0")->getData();

    EXPECT_EQ(testvalue1, testvalue2);

    testvalue1 =
      ds1->getRoot()->getGroup("fields2")->getGroup("b")->getView("i1")->getData();
    testvalue2 =
      ds2->getRoot()->getGroup("fields2")->getGroup("b")->getView("i1")->getData();

    EXPECT_EQ(testvalue1, testvalue2);

    delete ds1;
    delete ds2;
  }
}

//------------------------------------------------------------------------------
TEST(spio_serial, write_read_write)
{
  int my_rank, num_ranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  int num_files = std::max(num_ranks / 2, 1);
  std::stringstream sstr;
  sstr << "out_spio_WRW_" << num_ranks;
  std::string filename = sstr.str();

  // Initialize a datastore and dump to disk
  DataStore* ds = new DataStore();
  ds->getRoot()->createViewScalar("grp/i", 2);
  ds->getRoot()->createViewScalar("grp/f", 3.0);
  IOManager writer_a(MPI_COMM_WORLD);
  writer_a.write(ds->getRoot(), num_files, filename, PROTOCOL);

  // Create another DataStore to read into.
  DataStore ds_r;
  IOManager reader(MPI_COMM_WORLD);
  reader.read(ds_r.getRoot(), filename + ROOT_EXT);

  // Dump this datastore to disk.
  // Regression for sidre_hdf5 protocol:
  // This used to produce the following HDF5 error:
  //  HDF5-DIAG: Error detected in HDF5 (1.8.16) thread 0:
  //    #000: H5F.c line 522 in H5Fcreate(): unable to create file
  //      major: File accessibility
  //      minor: Unable to open file
  //    #001: H5Fint.c line 1024 in H5F_open(): unable to truncate a file which
  // is already open
  //      major: File accessibility
  //      minor: Unable to open file
  IOManager writer_b(MPI_COMM_WORLD);
  writer_b.write(ds_r.getRoot(), num_files, filename, PROTOCOL);
}
