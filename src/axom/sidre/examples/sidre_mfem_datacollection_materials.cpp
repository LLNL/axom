// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * @file sidre_mfem_datacollection_vis.cpp
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
  axom::sidre::MFEMSidreDataCollection dc("sidre_mfem_datacoll_materials_ex",
                                          &mesh);

  // _sidredc_material_matset_start
  // Inform the DataCollection that volume fraction information for a material set
  // called "matset" will be in fields called "volume_fraction"
  dc.AssociateMaterialSet("volume_fraction", "matset");
  // _sidredc_material_matset_end

  // _sidredc_material_depfield_start
  // Inform the DataCollection of a material-dependent field "density" associated
  // with a material set called "matset"
  dc.AssociateMaterialDependentField("density", "matset");
  // _sidredc_material_depfield_end

  // _sidredc_material_specset_start
  // Inform the DataCollection that species set data for the species set "specset"
  // associated with material set "matset" will be in fields called "partial_density"
  const bool volume_dependent = false;
  dc.AssociateSpeciesSet("partial_density", "specset", "matset", volume_dependent);
  // _sidredc_material_specset_end

  mfem::GridFunction vol_frac_1(&fes);
  vol_frac_1 = 0.65;
  mfem::GridFunction vol_frac_2(&fes);
  vol_frac_2 = 0.35;
  // _sidredc_material_matset_register_start
  dc.RegisterField("volume_fraction_001", &vol_frac_1);
  dc.RegisterField("volume_fraction_002", &vol_frac_2);
  // _sidredc_material_matset_register_end

  mfem::GridFunction density_indenpendent(&fes);
  density_indenpendent = 1.0;
  mfem::GridFunction density_dependent_1(&fes);
  density_dependent_1 = 0.2;
  mfem::GridFunction density_dependent_2(&fes);
  density_dependent_2 = 0.3;
  // _sidredc_material_depfield_register_start
  dc.RegisterField("density", &density_indenpendent);
  dc.RegisterField("density_001", &density_dependent_1);
  dc.RegisterField("density_002", &density_dependent_2);
  // _sidredc_material_depfield_register_end

  mfem::GridFunction partial_density_1_1(&fes);
  partial_density_1_1 = 0.4;
  mfem::GridFunction partial_density_1_2(&fes);
  partial_density_1_2 = 0.3;
  mfem::GridFunction partial_density_2_1(&fes);
  partial_density_2_1 = 0.2;
  mfem::GridFunction partial_density_2_2(&fes);
  partial_density_2_2 = 0.1;
  // _sidredc_material_specset_register_start
  dc.RegisterField("partial_density_001_001", &partial_density_1_1);
  dc.RegisterField("partial_density_001_002", &partial_density_1_2);
  dc.RegisterField("partial_density_002_001", &partial_density_2_1);
  dc.RegisterField("partial_density_002_002", &partial_density_2_2);
  // _sidredc_material_specset_register_end

  // This is where the time-dependent operator would be set up...

  // Save the initial state
  dc.SetCycle(0);   // Iteration counter
  dc.SetTime(0.0);  // Simulation time
  // Filename and protocol, both of which are optional
  dc.Save("sidre_mfem_datacoll_materials_ex", "sidre_conduit_json");
}
