// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**************************************************************************
 *************************************************************************/

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/sidre.hpp"

#include "mpi.h"
#include "CLI11/CLI11.hpp"

#ifdef AXOM_USE_SCR
  #include "scr.h"
#endif

using axom::sidre::DataStore;
using axom::sidre::DataType;
using axom::sidre::Group;
using axom::sidre::IOManager;
using namespace axom::utilities;

/** Define structure of datastore */
void initializeDS(int my_rank, DataStore* ds)
{
  Group* root = ds->getRoot();

  Group* flds = root->createGroup("fields");
  Group* ga = flds->createGroup("a");
  ga->createViewScalar<int>("i0", my_rank + 101);

  Group* flds2 = root->createGroup("fields2");
  Group* gb = flds2->createGroup("b");
  gb->createViewScalar<int>("i1", 4 * my_rank * my_rank + 404);

  return;
}

/** Write an output dataset via SCR */
bool dumpOutput(MPI_Comm comm,
                const std::string& file_base,
                int num_files,
                DataStore* ds)
{
  // Tell SCR we're starting an output dataset.
  // Each SCR output set should be given a name.
  // This name should be user-friendly, because
  // end users may need to type it at a command line.
  // It should also encode enough information that
  // one can construct the full path to the file on the
  // parallel file system.
  // SCR_Start_output must be called by all processes in MPI_COMM_WORLD.
  SCR_Start_output(file_base.c_str(), SCR_FLAG_OUTPUT);

  // write out the Datastore
  IOManager writer(comm, true);
  writer.write(ds->getRoot(), num_files, file_base, "sidre_hdf5");

  // Tell SCR this process has completed its checkpoint.
  // Each process should set valid=1 if it wrote its portion
  // of the checkpoint successfully and valid=0 otherwise.
  // SCR executes an allreduce across processes to identify
  // whether all ranks succeeded.
  // If any rank failed, SCR considers the checkpoint
  // to be invalid.
  // The return code will be SCR_SUCCESS if all procs
  // succeeded.
  // SCR_Complete_output must be called by all processes in MPI_COMM_WORLD.
  int valid = 1;
  int rc = SCR_Complete_output(valid);
  return (rc == SCR_SUCCESS);
}

/** Simple structure to hold the parsed command line arguments */
struct CommandLineArguments
{
  int m_numFiles;
  std::string m_fileBase;

  CommandLineArguments() : m_numFiles(0), m_fileBase("test.hdf") { }

  void parse(int argc, char** argv, CLI::App& app);
};

/** Parse the command line arguments */
void CommandLineArguments::parse(int argc, char** argv, CLI::App& app)
{
  app.add_option("-n,--num", m_numFiles, "Number of files per checkpoint")
    ->check(CLI::PositiveNumber);

  app.add_option("-f,--file", m_fileBase, "Base name of checkpoint files");

  app.get_formatter()->column_width(35);

  // Could throw an exception
  app.parse(argc, argv);
}

/** Terminates execution */
void quitProgram(int exitCode = 0)
{
  MPI_Finalize();
  exit(exitCode);
}

/**************************************************************************
 * Subroutine:  main
 * Purpose   :  This demonstrates use of SCR with SPIO to write an output
 *              datastore.
 *************************************************************************/

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  axom::slic::SimpleLogger logger;

  int my_rank, num_ranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  // parse the command line arguments
  CommandLineArguments args;
  CLI::App app {"SCR Output example"};

  try
  {
    args.parse(argc, argv, app);
  }
  catch(const CLI::ParseError& e)
  {
    int retval = -1;
    if(my_rank == 0)
    {
      retval = app.exit(e);
    }
    MPI_Bcast(&retval, 1, MPI_INT, 0, MPI_COMM_WORLD);
    quitProgram(retval);
  }

  // Extract values from command line.
  int num_files = args.m_numFiles;
  std::string file_base = args.m_fileBase;

  // Default to write a file per process.
  if(num_files == 0)
  {
    num_files = num_ranks;
  }

  // SCR_Init must be called by all processes in MPI_COMM_WORLD,
  // and it must be called after MPI_Init.
  SCR_Init();

  DataStore* ds = new DataStore();
  SLIC_ASSERT(ds);

  // Initialize DataStore
  initializeDS(my_rank, ds);

  // Write checkpoint via SCR.
  // One must write with MPI_COMM_WORLD when using SCR.
  // Also, the number of files should be the same as
  // the number of ranks in the job.
  dumpOutput(MPI_COMM_WORLD, file_base, num_files, ds);

  delete ds;

  // SCR_Finalize must be called by all processes in MPI_COMM_WORLD,
  // and it must be called before MPI_Finalize.
  // Among other things, SCR_Finalize flushes cached datasets to the
  // parallel file system.  It also signals the SCR run scripts
  // that the job should *not* be restarted within the current job
  // allocation.
  SCR_Finalize();

  MPI_Finalize();

  return 0;
}
