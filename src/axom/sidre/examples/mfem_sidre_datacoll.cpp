// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * @file mfem_sidre_datacoll.cpp
 * @brief This example code is a basic demonstration of Sidre's
 * MFEMSidreDataCollection class for visualizing a time-marching
 * simulation.  It is a more thorough version of the snippet
 * provided in Sidre's Sphinx documentation.
 */

// Datacollection header
#include "axom/sidre/core/MFEMSidreDataCollection.hpp"

// MFEM includes - needed to set up simulation
#include "mfem.hpp"

int main()
{
  // Built a 2D mesh with 100 square elements
  mfem::Mesh mesh(10, 10, mfem::Element::QUADRILATERAL);

  // Set up the FiniteElementSpace - needed for the grid functions
  // Initialize with H1 elements of order 1
  mfem::H1_FECollection fec(/*order=*/1, /*dimension=*/2);
  mfem::FiniteElementSpace fes(&mesh, &fec);

  // Initialize the datacollection with the mesh
  // Note: all fields (added with RegisterField) must be on this mesh
  axom::sidre::MFEMSidreDataCollection dc("mfem_sidre_datacoll_ex", &mesh);

  // This is where the time-dependent operator would be set up...

  // Initialize the solution field
  mfem::GridFunction soln(&fes);
  soln = 0.0;

  // Note: Any number of fields can be registered
  dc.RegisterField("solution", &soln);

  // Save the initial state
  dc.SetCycle(0);   // Iteration counter
  dc.SetTime(0.0);  // Simulation time
  // Filename and protocol, both of which are optional
  dc.Save("mfem_sidre_datacoll_ex", "sidre_hdf5");

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
    dc.Save("mfem_sidre_datacoll_ex", "sidre_hdf5");
  }
}
