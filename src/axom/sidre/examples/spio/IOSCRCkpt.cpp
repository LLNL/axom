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

  SLIC_ERROR_IF(argc != 3,
                "Missing command line arguments. \n\t"
                << "Usage: spio_IOCkpt <num_ckpts> <base_file_name>");

  int num_ckpts;
  std::string file_base;
  if (argc == 3)
  {
    num_ckpts = atoi(argv[1]);
    file_base = argv[2];
  }
  else
  {
    return 0;
  }

  int my_rank, num_ranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  int t_start = 0;

  // SCR restart loop
  // The application can continue trying to read a checkpoint
  // until it succeeds or it exhausts all available checkpoints.
  int have_restart = 0;
  int restarted = 0;
  do {
    // Ask SCR if there is a checkpoint to read.
    // SCR loads the next most recent checkpoint that
    // it has not internally marked as "failed".
    // If there is no checkpoint, SCR will set have_restart to 0.
    // SCR will set have_restart to 1 if a checkpoint is
    // available, and copy the name of that checkpoint into ckptname.
    // This name is the same string the application provided when
    // it created the checkpoint in the work loop below.
    char ckptname[SCR_MAX_FILENAME];
    SCR_Have_restart(&have_restart, ckptname);
    if (have_restart) {
      // Tell SCR we're starting to read the checkpoint.
      // One must only call SCR_Start_restart if SCR_Have_restart
      // reports that there is a checkpoint to restart from.
      // SCR returns the same name in ckptname as SCR_Have_restart.
      SCR_Start_restart(ckptname);

      // It is common to use the checkpoint name to compute
      // the path to the checkpoint file.
      // One should build the path to the file on the parallel
      // file system.
      std::string root_file = std::string(ckptname) + "." + file_base + ".root";

      DataStore* ds = new DataStore();
      SLIC_ASSERT(ds);
      Group* root = ds->getRoot();
    
      IOManager reader(MPI_COMM_WORLD, true);
      reader.read(root, root_file, false, true);
    
      int timestep = root->getView("fields/a/timestep")->getScalar();

      delete ds;

      // Tell SCR whether this process succeeded reading its checkpoint.
      // Each process should set valid=1 if it read its portion
      // of the checkpoint successfully and valid=0 otherwise.
      // SCR executes an allreduce across processes to identify
      // whether all ranks succeeded.
      // If any rank failed, SCR considers the checkpoint
      // to be invalid, and then SCR will load the next most
      // recent checkpoint if one exists.
      // The return code will be SCR_SUCCESS if all procs
      // succeeded.
      int valid = 1;
      int rc = SCR_Complete_restart(valid);
      restarted = (rc == SCR_SUCCESS);

      // If we succeeded, then update initial state
      if (restarted) {
        t_start = timestep + 1;
      }
    }
  } while (have_restart && !restarted);

  // Application work loop
  int t;
  for (t = t_start; t < t_start + num_ckpts; t++) {
    /////////////////////////
    // Do actual work ...
    /////////////////////////

    // Ask SCR whether it's time to checkpoint is optional.
    // For checkpoints that are purely defensive, SCR can
    // help guide the application to the proper checkpoint
    // frequency for exampled based on checkpoint cost and
    // expected failure frequency.
    int need_checkpoint = 0;
    SCR_Need_checkpoint(&need_checkpoint);
    if (need_checkpoint) {

      // Tell SCR we're starting a checkpoint.
      // Each SCR output set should be given a name.
      // This name should be user-friendly, because
      // end users may need to type it at a command line.
      // It should also encode enough information that
      // one can construct the full path to the file on the
      // parallel file system.
      std::string ckptname = "time." + std::to_string(t);
      SCR_Start_output(ckptname.c_str(), SCR_FLAG_CHECKPOINT);

      DataStore* ds = new DataStore();
      SLIC_ASSERT(ds);
      Group* root = ds->getRoot();

      Group* flds = root->createGroup("fields");
      Group* flds2 = root->createGroup("fields2");

      Group* ga = flds->createGroup("a");
      Group* gb = flds2->createGroup("b");
      ga->createViewScalar<int>("timestep", t);
      ga->createViewScalar<int>("i0", my_rank + 101);
      gb->createViewScalar<int>("i1", 4*my_rank*my_rank + 404);

      std::string ckpt_path = ckptname + "." + file_base;

      // One must write with MPI_COMM_WORLD when using SCR.
      // Also, the number of files should be the same as
      // the number of ranks in the job.
      IOManager writer(MPI_COMM_WORLD, true);
      writer.write(root, num_ranks, ckpt_path, "sidre_hdf5");

      MPI_Barrier(MPI_COMM_WORLD);

      delete ds;

      // Tell SCR this process has completed its checkpoint.
      // Each process should set valid=1 if it wrote its portion
      // of the checkpoint successfully and valid=0 otherwise.
      // SCR executes an allreduce across processes to identify
      // whether all ranks succeeded.
      // If any rank failed, SCR considers the checkpoint
      // to be invalid.
      // The return code will be SCR_SUCCESS if all procs
      // succeeded.
      int valid = 1;
      int complete_rc = SCR_Complete_output(valid);
      if (complete_rc != SCR_SUCCESS) {
        // some process failed to checkpoint
      }

      // An application can optionally ask SCR whether it should
      // stop.  This is important when datasets are cached in order
      // to leave SCR time to flush those datasets to the parallel file
      // system before the job's allocatino expires.
      int should_exit = 0;
      SCR_Should_exit(&should_exit);
      if (should_exit)
      {
          break;
      }
    } // need checkpoint
  }

  SCR_Finalize();
  MPI_Finalize();

  return 0;
}
