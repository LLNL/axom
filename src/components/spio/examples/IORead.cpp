/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */


/**************************************************************************
 *************************************************************************/

#include "mpi.h"

#include "sidre/DataGroup.hpp"
#include "sidre/DataStore.hpp"
#include "spio/IOManager.hpp"

using asctoolkit::sidre::DataGroup;
using asctoolkit::sidre::DataStore;
using asctoolkit::spio::IOManager;

/**************************************************************************
 * Subroutine:  main
 * Purpose   :
 *************************************************************************/

int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);

  SLIC_ASSERT(argc == 2);

  DataStore * ds = new DataStore();
  SLIC_ASSERT(ds);
  DataGroup * root = ds->getRoot();

  size_t num_files = 0;
  std::string root_file;
  if (argc == 2) {
    root_file = argv[1];
    hid_t root_file_id = H5Fopen(root_file.c_str(),
                                 H5F_ACC_RDWR,
                                 H5P_DEFAULT);
    SLIC_ASSERT(root_file_id >= 0);

    hid_t num_files_id = H5Dopen(root_file_id, "num_files", H5P_DEFAULT);
    SLIC_ASSERT(num_files_id >= 0);

    int tmp_num_files;
    herr_t errv = H5Dread(num_files_id,
            H5T_NATIVE_INT,
            H5S_ALL,
            H5S_ALL,
            H5P_DEFAULT,
            &tmp_num_files);
    SLIC_ASSERT(errv >= 0);
    SLIC_ASSERT(tmp_num_files >= 0);
    num_files = static_cast<size_t>(tmp_num_files);
  } else {
    return 0;
  }

  std::vector<DataGroup *> groups;
  groups.push_back(root);

  IOManager reader(MPI_COMM_WORLD, &(groups[0]), groups.size(), num_files);
  reader.read(root_file);

  delete ds;

  MPI_Finalize();


  return 0;
}
