// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
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

// MFEM includes - needed to set up simulation
#include "mfem.hpp"

#ifdef AXOM_USE_MPI
  #include "mpi.h"
#endif

int main(int argc, char* argv[])
{
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  MPI_Init(&argc, &argv);
#endif

  mfem::Mesh* mesh = nullptr;

  // Built a 2D mesh with 100 square elements
  mfem::Mesh serial_mesh(10, 10, mfem::Element::QUADRILATERAL);
  mesh = &serial_mesh;

#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  mfem::ParMesh parallel_mesh(MPI_COMM_WORLD, serial_mesh);
  mesh = &parallel_mesh;
#endif

  // Set up the FiniteElementSpace - needed for the grid functions
  // Initialize with H1 elements of order 1
  mfem::H1_FECollection fec(/*order=*/1, /*dimension=*/2);
  mfem::FiniteElementSpace fes(mesh, &fec);

  // Initialize the datacollection with the mesh
  // Note: all fields (added with RegisterField) must be on this mesh
  axom::sidre::MFEMSidreDataCollection dc("sidre_mfem_datacoll_ex", mesh);

// This is where the time-dependent operator would be set up...

// Initialize the solution field
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  mfem::ParFiniteElementSpace par_fes(fes, parallel_mesh);
  mfem::ParGridFunction soln(&par_fes);
#else
  mfem::GridFunction soln(&fes);
#endif
  soln = 0.0;

  // Note: Any number of fields can be registered
  dc.RegisterField("solution", &soln);

  // Save the initial state
  dc.SetCycle(0);   // Iteration counter
  dc.SetTime(0.0);  // Simulation time
  // Filename and protocol, both of which are optional
  dc.Save("sidre_mfem_datacoll_ex", "sidre_hdf5");

  // Sample time parameters
  int n_iter = 10;
  double dt = 0.05;

  for(int i = 0; i < n_iter; i++)
  {
    // Calculate the next iteration of the solution field...
    // For simplicity, every element in the field is set to the current time
    soln = dt * i;

    // then save it after updating the time information...
    dc.SetCycle(i);
    dc.SetTime(dt * i);
    dc.Save("sidre_mfem_datacoll_ex", "sidre_hdf5");
  }

#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  MPI_Finalize();
#endif
}
