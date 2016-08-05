/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */


/*!
 * \file
 *
 * Simple reproducer for bug in Spio dealing with consecutive
 * reading and writing on same set of files.
 */

#include "spio/IOManager.hpp"
#include "sidre/sidre.hpp"

#include "mpi.h"

#include <string>

using axom::spio::IOManager;
using axom::sidre::Group;
using axom::sidre::DataStore;
using axom::sidre::DataType;
using axom::sidre::View;

const std::string PROTOCOL = "sidre_hdf5";
const std::string ROOT_EXT = ".root";


/*!
 * Dump a Sidre Group using Spio
 */
void dump(const std::string& fName, Group* grp, int numFiles = 1)
{
    IOManager writer(MPI_COMM_WORLD);
    writer.write(grp, numFiles, fName, PROTOCOL);
}

/*!
 * Load a Sidre Group using Spio
 */
void load(const std::string& fName, Group* grp)
{
    IOManager reader(MPI_COMM_WORLD);
    reader.read(grp, fName + ROOT_EXT);
}


//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);

  int my_rank, num_ranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  int num_files = std::max( num_ranks / 2, 1);
  std::string filename = "out_spio_RWR";


  DataStore * ds = new DataStore();
  ds->getRoot()->createViewScalar("grp/i",2);
  ds->getRoot()->createViewScalar("grp/f",3.0);

  dump(filename, ds->getRoot(), num_files);

  // Create another DataStore to read into.
  DataStore ds_r;

  load(filename, ds_r.getRoot());

  // BUG HERE.  The following dump() produces an HDF5 error
  //    When loading the root file
  //  HDF5-DIAG: Error detected in HDF5 (1.8.16) thread 0:
  //    #000: H5F.c line 522 in H5Fcreate(): unable to create file
  //      major: File accessibility
  //      minor: Unable to open file
  //    #001: H5Fint.c line 1024 in H5F_open(): unable to truncate a file which is already open
  //      major: File accessibility
  //      minor: Unable to open file
  dump(filename, ds_r.getRoot(), num_files);


  MPI_Finalize();

  return 0;
}

