/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
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

/**************************************************************************
 *************************************************************************/

#include "mpi.h"

#include "slic/slic.hpp"
#include "slic/UnitTestLogger.hpp"

#include "sidre/Group.hpp"
#include "sidre/DataStore.hpp"
#include "sidre/IOManager.hpp"

using axom::sidre::Group;
using axom::sidre::DataStore;
using axom::sidre::IOManager;

/**************************************************************************
 * Subroutine:  main
 * Purpose   :
 *************************************************************************/

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  axom::slic::UnitTestLogger logger;

  SLIC_ERROR_IF(argc != 2,
                "Missing required command line argument. \n\t"
                << "Usage: spio_IORead <sidre_root_file>");

  DataStore* ds = new DataStore();
  SLIC_ASSERT(ds);
  Group* root = ds->getRoot();

  std::string root_file;
  if (argc == 2)
  {
    root_file = argv[1];
  }
  else
  {
    return 0;
  }

  IOManager reader(MPI_COMM_WORLD);
  reader.read(root, root_file);

  delete ds;

  MPI_Finalize();


  return 0;
}
