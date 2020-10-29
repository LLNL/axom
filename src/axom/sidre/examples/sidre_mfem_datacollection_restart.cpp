// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * @file sidre_mfem_datacollection.cpp
 * @brief This example code is a basic demonstration of Sidre's
 * MFEMSidreDataCollection class for visualizing a time-marching
 * simulation.  It is a more thorough version of the snippet
 * provided in Sidre's Sphinx documentation.
 */

// Datacollection header
#include "axom/sidre/core/MFEMSidreDataCollection.hpp"

#include <memory>

#ifdef AXOM_USE_MPI
  #include "mpi.h"
#endif

// MFEM includes - needed to set up simulation
#include "mfem.hpp"

#include "axom/core/utilities/FileUtilities.hpp"

// Abstraction for an object that may or may not retain ownership of an object
// via an owning- or non-owning-pointer
template <typename T>
class MaybeOwner
{
public:
  // Would also be useful to have a make_MaybeOwner ctor that forwards args
  MaybeOwner& operator=(std::unique_ptr<T>&& ptr)
  {
    m_non_owning = nullptr;
    m_owning = std::move(ptr);
    return *this;
  }
  // Call move constructor - really doesn't need to be on the heap though
  // Use std::optional when available
  MaybeOwner& operator=(T&& obj)
  {
    m_non_owning = nullptr;
    m_owning.reset(new T(std::move(obj)));  // Calls move constructor
    return *this;
  }
  // Avoid null check
  MaybeOwner& operator=(T& obj)
  {
    m_owning.reset();
    m_non_owning = &obj;
    return *this;
  }

  MaybeOwner& operator=(T* ptr)
  {
    SLIC_ERROR_IF(!ptr, "Attempted to set a MaybeOwner to a nullptr");
    m_owning.reset();
    m_non_owning = ptr;
    return *this;
  }

  T& get()
  {
    checkState();
    return (m_owning) ? *m_owning : *m_non_owning;
  }

  const T& get() const
  {
    checkState();
    return (m_owning) ? *m_owning : *m_non_owning;
  }

  operator T&() { return get(); }
  operator const T&() const { return get(); }

private:
  std::unique_ptr<T> m_owning;  // "nullptr" - value-initialized ptr member - by default
  T* m_non_owning = nullptr;
  void checkState()
  {
    SLIC_ERROR_IF(m_owning && m_non_owning,
                  "Invalid MaybeOwner: Both pointers exist");
    SLIC_ERROR_IF(!m_owning && !m_non_owning,
                  "Invalid MaybeOwner: Neither pointer exists");
  }
};

// Stores the state of the simulation - a mesh, fields, and associated objects
struct SimulationState
{
  SimulationState(axom::sidre::MFEMSidreDataCollection& dc) : datacoll(dc)
  {
    // FIXME: Hardcoded, may not even work for a restart of a restart
    const std::string saved_data_filename =
      dc.GetCollectionName() + "_000000.root";

    // Check if this is a restart run
    if(axom::utilities::filesystem::pathExists(saved_data_filename))
    {
      SLIC_INFO("Reading in existing data and restarting...");
      // If it is, we can load everything in and "unwrap" to fill in the state
      dc.Load();
      // The Mesh, GridFunction, etc, objects already exist and can be accessed
      mesh = dc.GetMesh();
      soln_field = dc.GetField("solution");
      fespace = soln_field.get().FESpace();
      fecoll = fespace.get().FEColl();
    }
    // Otherwise it's a nominal run so we have to create everything
    // In a realistic simulation this is where an input file might be used
    else
    {
      SLIC_INFO("Starting a new simulation");
      // Everything gets created on the heap so its lifetime can be managed by
      // the MaybeOwners - needs to persist after this function

      // Built a 2D mesh with 100 square elements
      std::unique_ptr<mfem::Mesh> tmp_mesh(
        new mfem::Mesh(10, 10, mfem::Element::QUADRILATERAL));

#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
      tmp_mesh.reset(new mfem::ParMesh(MPI_COMM_WORLD, *tmp_mesh));
#endif
      mesh = std::move(tmp_mesh);
      // Set up the DataCollection with the newly created mesh
      dc.SetMesh(&mesh.get());

      // Set up the FiniteElementSpace - needed for the grid functions
      // Initialize with H1 elements of order 1
      fecoll = std::unique_ptr<mfem::H1_FECollection>(
        new mfem::H1_FECollection(/*order=*/1, mesh.get().Dimension()));
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
      auto par_mesh = dynamic_cast<mfem::ParMesh*>(&mesh.get());
      fespace = std::unique_ptr<mfem::ParFiniteElementSpace>(
        new mfem::ParFiniteElementSpace(par_mesh, &fecoll.get()));
#else
      fespace = std::unique_ptr<mfem::FiniteElementSpace>(
        new mfem::FiniteElementSpace(&mesh.get(), &fecoll.get()));
#endif

      // Initialize the solution field
      std::unique_ptr<mfem::GridFunction> tmp_soln_field;

      // Set the data to nullptr so the datacollection will initialize it with
      // its own managed data (needed for a restart)
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
      auto par_fespace =
        dynamic_cast<mfem::ParFiniteElementSpace*>(&fespace.get());
      soln_field = std::unique_ptr<mfem::ParGridFunction>(
        new mfem::ParGridFunction(par_fespace, static_cast<double*>(nullptr)));
#else
      soln_field = std::unique_ptr<mfem::GridFunction>(
        new mfem::GridFunction(&fespace.get(), nullptr));
#endif
      dc.RegisterField("solution", &soln_field.get());

      // Intialize to zero as our "initial conditions"
      soln_field.get() = 0.0;

      // Set t = 0 state info
      dc.SetCycle(0);   // Iteration counter
      dc.SetTime(0.0);  // Simulation time
    }
  }

  void step(double dt)
  {
    // Update simulation state variables
    double t = datacoll.GetTime();
    t += dt;
    datacoll.SetTime(t);

    const int cycle = datacoll.GetCycle();
    datacoll.SetTime(cycle + 1);

    // Calculate the next iteration of the solution field...
    // For simplicity, every element in the field is set to the current time
    soln_field.get() = t;
  }

  MaybeOwner<mfem::Mesh> mesh;
  MaybeOwner<const mfem::FiniteElementCollection> fecoll;
  MaybeOwner<mfem::FiniteElementSpace> fespace;
  MaybeOwner<mfem::GridFunction> soln_field;
  axom::sidre::MFEMSidreDataCollection& datacoll;
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

  // Initialize the simulation data structures
  SimulationState sim_state(dc);

  // This is where the time-dependent operator would be set up...

  // Save initial state of simulation
  dc.Save("sidre_mfem_datacoll_restart_ex", "sidre_hdf5");

  // Sample time parameters
  int n_iter = 10;
  double dt = 0.05;

  for(int i = 0; i < n_iter; i++)
  {
    sim_state.step(dt);

    // then save it
    dc.Save("sidre_mfem_datacoll_restart_ex", "sidre_hdf5");
  }

#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  MPI_Finalize();
#endif
}
