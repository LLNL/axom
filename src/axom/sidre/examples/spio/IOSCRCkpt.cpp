// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**************************************************************************
 *************************************************************************/

#include "mpi.h"

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/sidre.hpp"

#ifdef AXOM_USE_SCR
#include "scr.h"
#endif

using axom::sidre::Group;
using axom::sidre::DataStore;
using axom::sidre::DataType;
using axom::sidre::IOManager;
using namespace axom::utilities;

/**************************************************************************
 * Subroutine:  main
 * Purpose   :
 *************************************************************************/

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  SCR_Init();
  axom::slic::UnitTestLogger logger;

  SLIC_ERROR_IF(argc != 2,
                "Missing command line arguments. \n\t"
                << "Usage: spio_IOWrite <base_file_name>");

  std::string file_base;
  if (argc == 2)
  {
    file_base = argv[1];
  }
  else
  {
    return 0;
  }

  int my_rank, num_ranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  // restart loop
  int have_restart = 0;
  int restarted = 0;
  do {
    // ask SCR if there is a checkpoint to read
    char ckptname[SCR_MAX_FILENAME];
    SCR_Have_restart(&have_restart, ckptname);
    if (have_restart) {
      // tell SCR we're starting to read the checkpoint
      SCR_Start_restart(ckptname);

      //std::string root_file = std::string(ckptname) + "." + file_base;
      std::string root_file = std::string(ckptname) + "." + file_base + ".root";

      DataStore* ds = new DataStore();
      SLIC_ASSERT(ds);
      Group* root = ds->getRoot();
    
      IOManager reader(MPI_COMM_WORLD, true);
      reader.read(root, root_file, false, true);
    
      delete ds;

      // tell SCR whether this process succeeded,
      // return code tells us whether all procs succeeded
      int rc = SCR_Complete_restart(1);
      restarted = (rc == SCR_SUCCESS);
    }
  } while (have_restart && !restarted);

  int t;
  for (t = 0; t < 10; t++) {
    // call to ask whether it's time to checkpoint is optional
    int need_checkpoint = 0;
    SCR_Need_checkpoint(&need_checkpoint);
    //if (need_checkpoint) {

      // tell SCR we're starting a checkpoint
      std::string ckptname = "time." + std::to_string(t);
      SCR_Start_output(ckptname.c_str(), SCR_FLAG_CHECKPOINT);

      DataStore* ds = new DataStore();
      SLIC_ASSERT(ds);
      Group* root = ds->getRoot();

      Group* flds = root->createGroup("fields");
      Group* flds2 = root->createGroup("fields2");

      Group* ga = flds->createGroup("a");
      Group* gb = flds2->createGroup("b");
      //ga->createViewScalar<int>("timestep", timestep);
      ga->createViewScalar<int>("i0", my_rank + 101);
      gb->createViewScalar<int>("i1", 4*my_rank*my_rank + 404);

#if 0
      if (my_rank == 0)
      {
        std::string dir;
        filesystem::getDirName(dir, file_base);
        filesystem::makeDirsForPath(dir);
      }
      MPI_Barrier(MPI_COMM_WORLD);
#endif

      std::string ckpt_path = ckptname + "." + file_base;

      IOManager writer(MPI_COMM_WORLD, true);
      writer.write(root, num_ranks, ckpt_path, "sidre_hdf5");

      MPI_Barrier(MPI_COMM_WORLD);

      delete ds;

      // tell SCR this process has completed its checkpoint
      // return code tells us whether all procs succeeded
      SCR_Complete_output(1);

      // call to ask whether it's time to exit is optional
      int should_exit = 0;
      SCR_Should_exit(&should_exit);
      //if (should_exit) break;
    //} // need checkpoint
  }

  SCR_Finalize();
  MPI_Finalize();

  return 0;
}
