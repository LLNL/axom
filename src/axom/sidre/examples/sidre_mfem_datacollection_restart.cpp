// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * @file sidre_mfem_datacollection_restart.cpp
 * @brief This example code is a basic demonstration of Sidre's
 * MFEMSidreDataCollection class for restarting a simulation.
 * 
 * To reload from a file, the cycle number to reload can be
 * specified as a command-line argument.
 * 
 * For example, run
 * @code{.sh}
 * ./sidre_mfem_datacollection_restart # generates cycles 0-9
 * # then...
 * ./sidre_mfem_datacollection_restart 9 # loads cycle 9 and continues
 * @endcode
 */

// Datacollection header
#include "axom/sidre/core/MFEMSidreDataCollection.hpp"

#include <memory>  // for unique_ptr

#ifdef AXOM_USE_MPI
  #include "mpi.h"
#endif

// MFEM includes - needed to set up simulation
#include "mfem.hpp"
#include "CLI11/CLI11.hpp"

// Stores the state of the simulation - a mesh, fields, and associated objects
class SimulationState
{
public:
  // Initializes the simulation, using an existing file with specified cycle if
  // cycle_to_load is specified (>= 0)
  SimulationState(axom::sidre::MFEMSidreDataCollection& dc,
                  const int cycle_to_load)
    : datacoll(dc)
  {
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
    MPI_Comm_rank(datacoll.GetComm(), &m_rank);
#endif
    // Check if this is a restart run
    if(cycle_to_load >= 0)
    {
      // If it is, we can load everything in and "unwrap" to fill in the state
      datacoll.Load(cycle_to_load);
      reloadSim();
    }
    // Otherwise it's a nominal run so we have to create everything
    // In a realistic simulation this is where an input file might be used
    else
    {
      setupNewSim();
    }
  }

  ~SimulationState()
  {
    if(owns_data)
    {
      delete mesh;
      delete fecoll;
      delete fespace;
      delete soln_field;
    }
  }

  // A simulated "step" of the simulation
  void step(double dt)
  {
    // Update simulation state variables
    double t = datacoll.GetTime();
    t += dt;
    datacoll.SetTime(t);

    const int cycle = datacoll.GetCycle();
    datacoll.SetCycle(cycle + 1);

    // Calculate the next iteration of the solution field...
    // For simplicity, every element in the field is set to the current time
    *soln_field = t;
  }

private:
  // Simulation state setup
  void setupNewSim()
  {
    SLIC_INFO_IF(!m_rank, "Starting a new simulation");
    // Everything here is managed by the SimulationState object
    owns_data = true;

    // Build a 2D mesh with 100 square elements
    mesh = new mfem::Mesh(10, 10, mfem::Element::QUADRILATERAL);

#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
    mfem::Mesh* tmp_mesh = mesh;
    mesh = new mfem::ParMesh(MPI_COMM_WORLD, *tmp_mesh);
    delete tmp_mesh;
#endif
    // Set up the DataCollection with the newly created mesh
    datacoll.SetMesh(mesh);

    // Set up the FiniteElementSpace - needed for the grid functions
    // Initialize with H1 elements of order 1
    fecoll = new mfem::H1_FECollection(/*order=*/1, mesh->Dimension());
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
    auto par_mesh = dynamic_cast<mfem::ParMesh*>(mesh);
    fespace = new mfem::ParFiniteElementSpace(par_mesh, fecoll);
#else
    fespace = new mfem::FiniteElementSpace(mesh, fecoll);
#endif

    // Initialize the solution field

    // Set the data to nullptr so the datacollection will initialize it with
    // its own managed data (needed for a restart)
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
    auto par_fespace = dynamic_cast<mfem::ParFiniteElementSpace*>(fespace);
    soln_field =
      new mfem::ParGridFunction(par_fespace, static_cast<double*>(nullptr));
#else
    soln_field = new mfem::GridFunction(fespace, nullptr);
#endif
    datacoll.RegisterField("solution", soln_field);

    // Intialize to zero as our "initial conditions"
    *soln_field = 0.0;

    // Set t = 0 state info
    datacoll.SetCycle(0);   // Iteration counter
    datacoll.SetTime(0.0);  // Simulation time
  }

  // Sets up the state with non-owning pointers
  void reloadSim()
  {
    SLIC_INFO_IF(!m_rank,
                 "Reading in existing data and restarting from iteration "
                   << datacoll.GetCycle() << " at time " << datacoll.GetTime());
    // The Mesh, GridFunction, etc, objects already exist and can be accessed
    mesh = datacoll.GetMesh();
    soln_field = datacoll.GetField("solution");
    fespace = soln_field->FESpace();
    fecoll = fespace->FEColl();
  }

  // FEM-related objects needed as part of a simulation
  // In a real simulation these would be exposed via accessors
  mfem::Mesh* mesh;
  const mfem::FiniteElementCollection* fecoll;
  mfem::FiniteElementSpace* fespace;
  mfem::GridFunction* soln_field;
  bool owns_data = false;

  // A reference to the datacollection so it can be updated with time/cycle
  // information on each time step
  axom::sidre::MFEMSidreDataCollection& datacoll;

  // The MPI rank, used to display messages on just one node
  int m_rank = 0;
};

int main(int argc, char* argv[])
{
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  MPI_Init(&argc, &argv);
#endif

  // Initialize the datacollection
  // Needs to be configured to own the mesh data so all mesh data is saved to datastore/output file
  const bool owns_mesh_data = true;
  axom::sidre::MFEMSidreDataCollection dc("sidre_mfem_datacoll_restart_ex",
                                          nullptr,
                                          owns_mesh_data);
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  dc.SetComm(MPI_COMM_WORLD);
#endif

  // Command-line argument to load in a specific cycle - optional
  CLI::App app {"Example of Axom's MFEMSidreDataCollection for restarts"};
  int cycle_to_load = -1;
  app.add_option("--cycle", cycle_to_load, "Optional simulation cycle to load");
  CLI11_PARSE(app, argc, argv);

  // Initialize the simulation data structures
  SimulationState sim_state(dc, cycle_to_load);

  // This is where the time-dependent operator would be set up...

  // Save initial state of simulation
  dc.Save("sidre_mfem_datacoll_restart_ex", "sidre_hdf5");

  // Sample time parameters
  int n_iter = 10;
  double dt = 0.05;

  for(int i = 0; i < n_iter; i++)
  {
    sim_state.step(dt);
  }

  // then save it at the end of the simulation
  dc.Save("sidre_mfem_datacoll_restart_ex", "sidre_hdf5");

#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  MPI_Finalize();
#endif
}
