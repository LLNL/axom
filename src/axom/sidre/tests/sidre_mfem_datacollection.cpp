// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#ifndef AXOM_USE_MFEM
  #error This file requires MFEM
#endif

#include "mfem.hpp"

#include "gtest/gtest.h"

#include "axom/sidre/core/sidre.hpp"
#include "axom/sidre/core/MFEMSidreDataCollection.hpp"

using axom::sidre::Group;
using axom::sidre::MFEMSidreDataCollection;

const double EPSILON = 1.0e-6;

std::string testName()
{
  return ::testing::UnitTest::GetInstance()->current_test_info()->name();
}

TEST(sidre_datacollection, dc_alloc_no_mesh)
{
  MFEMSidreDataCollection sdc(testName());
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_alloc_owning_mesh)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  bool owns_mesh = true;
  MFEMSidreDataCollection sdc(testName(), &mesh, owns_mesh);
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_alloc_nonowning_mesh)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  bool owns_mesh = false;
  MFEMSidreDataCollection sdc(testName(), &mesh, owns_mesh);
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_register_empty_field)
{
  MFEMSidreDataCollection sdc(testName());
  mfem::GridFunction gf;
  sdc.RegisterField("test_field", &gf);
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_register_partial_field)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  mfem::H1_FECollection fec(2);
  mfem::FiniteElementSpace fes(&mesh, &fec);
  mfem::GridFunction gf(&fes);

  MFEMSidreDataCollection sdc(testName(), &mesh);
  sdc.RegisterField("test_field", &gf);
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_update_state)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  MFEMSidreDataCollection sdc(testName(), &mesh);

  // Arbitrary values for the "state" part of blueprint
  sdc.SetCycle(3);
  sdc.SetTime(1.1);
  sdc.SetTimeStep(0.4);
  // Force refresh from the superclass
  sdc.UpdateStateToDS();

  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_save)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  MFEMSidreDataCollection sdc(testName(), &mesh);

  sdc.Save();

  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_reload_gf)
{
  const std::string field_name = "test_field";
  // 2D mesh divided into triangles
  mfem::Mesh mesh(10, 10, mfem::Element::TRIANGLE);
  mfem::H1_FECollection fec(1, mesh.Dimension());
  mfem::FiniteElementSpace fes(&mesh, &fec);

  // The mesh and field(s) must be owned by Sidre to properly manage data in case of
  // a simulated restart (save -> load)
  bool owns_mesh = true;
  MFEMSidreDataCollection sdc_writer(testName(), &mesh, owns_mesh);
  mfem::GridFunction gf_write(&fes, nullptr);

  // Register to allocate storage internally, then write to it
  sdc_writer.RegisterField(field_name, &gf_write);

  mfem::ConstantCoefficient three_and_a_half(3.5);
  gf_write.ProjectCoefficient(three_and_a_half);

  EXPECT_TRUE(sdc_writer.verifyMeshBlueprint());

#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  sdc_writer.SetComm(MPI_COMM_WORLD);
#endif

  sdc_writer.SetCycle(0);
  sdc_writer.Save();

  // No mesh is used here
  MFEMSidreDataCollection sdc_reader(testName());

#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  sdc_reader.SetComm(MPI_COMM_WORLD);
#endif

  sdc_reader.Load();

  // No need to reregister, it already exists
  auto gf_read = sdc_reader.GetField(field_name);

  // Make sure the gridfunction was actually read in
  EXPECT_LT(gf_read->ComputeL2Error(three_and_a_half), EPSILON);

  EXPECT_TRUE(sdc_reader.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_reload_gf_vdim)
{
  const std::string field_name = "test_field";
  const int vdim = 2;
  // 2D mesh divided into triangles
  mfem::Mesh mesh(10, 10, mfem::Element::TRIANGLE);
  mfem::H1_FECollection fec(1, mesh.Dimension());
  mfem::FiniteElementSpace fes(&mesh, &fec, vdim, mfem::Ordering::byVDIM);

  // The mesh and field(s) must be owned by Sidre to properly manage data in case of
  // a simulated restart (save -> load)
  bool owns_mesh = true;
  MFEMSidreDataCollection sdc_writer(testName(), &mesh, owns_mesh);
  mfem::GridFunction gf_write(&fes, nullptr);

  // Register to allocate storage internally, then write to it
  sdc_writer.RegisterField(field_name, &gf_write);

  mfem::ConstantCoefficient three_and_a_half(3.5);
  gf_write.ProjectCoefficient(three_and_a_half);

  EXPECT_TRUE(sdc_writer.verifyMeshBlueprint());

#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  sdc_writer.SetComm(MPI_COMM_WORLD);
#endif

  sdc_writer.SetCycle(0);
  sdc_writer.Save();

  // No mesh is used here
  MFEMSidreDataCollection sdc_reader(testName());

#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  sdc_reader.SetComm(MPI_COMM_WORLD);
#endif

  sdc_reader.Load();

  // No need to reregister, it already exists
  auto gf_read = sdc_reader.GetField(field_name);

  // Make sure the gridfunction was actually read in
  EXPECT_LT(gf_read->ComputeL2Error(three_and_a_half), EPSILON);
  EXPECT_EQ(gf_read->FESpace()->GetOrdering(), mfem::Ordering::byVDIM);
  EXPECT_EQ(gf_read->FESpace()->GetVDim(), vdim);

  EXPECT_TRUE(sdc_reader.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_reload_mesh)
{
  const std::string field_name = "test_field";
  // 2D mesh divided into triangles
  mfem::Mesh mesh(10, 10, mfem::Element::TRIANGLE);
  mfem::H1_FECollection fec(1, mesh.Dimension());
  mfem::FiniteElementSpace fes(&mesh, &fec);

  // The mesh and field(s) must be owned by Sidre to properly manage data in case of
  // a simulated restart (save -> load)
  bool owns_mesh = true;
  MFEMSidreDataCollection sdc_writer(testName(), &mesh, owns_mesh);

  EXPECT_TRUE(sdc_writer.verifyMeshBlueprint());

  // Save some basic info about the mesh
  const int n_verts = sdc_writer.GetMesh()->GetNV();
  const int n_ele = sdc_writer.GetMesh()->GetNE();
  const int n_bdr_ele = sdc_writer.GetMesh()->GetNBE();

#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  sdc_writer.SetComm(MPI_COMM_WORLD);
#endif

  sdc_writer.SetCycle(0);
  sdc_writer.Save();

  // No mesh is used here
  MFEMSidreDataCollection sdc_reader(testName());

#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  sdc_reader.SetComm(MPI_COMM_WORLD);
#endif

  sdc_reader.Load();

  // Make sure the mesh was actually reconstructed
  EXPECT_EQ(sdc_reader.GetMesh()->GetNV(), n_verts);
  EXPECT_EQ(sdc_reader.GetMesh()->GetNE(), n_ele);
  EXPECT_EQ(sdc_reader.GetMesh()->GetNBE(), n_bdr_ele);

  EXPECT_TRUE(sdc_reader.verifyMeshBlueprint());
}

#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)

TEST(sidre_datacollection, dc_alloc_owning_parmesh)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  mfem::ParMesh parmesh(MPI_COMM_WORLD, mesh);
  bool owns_mesh = true;
  MFEMSidreDataCollection sdc(testName(), &parmesh, owns_mesh);
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_alloc_nonowning_parmesh)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  mfem::ParMesh parmesh(MPI_COMM_WORLD, mesh);
  bool owns_mesh = false;
  MFEMSidreDataCollection sdc(testName(), &parmesh, owns_mesh);
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

void dumpDataCollection(MFEMSidreDataCollection& sdc)
{
  conduit::Node mesh_node;
  sdc.GetBPGroup()->createNativeLayout(mesh_node);
  SLIC_INFO(mesh_node.to_yaml());
}

// Helper function to check for the existence of two sidre::Views and to check
// that they each refer to the same block of data
void checkReferentialEquality(axom::sidre::Group* grp,
                              const std::string& first,
                              const std::string& second)
{
  ASSERT_TRUE(grp->hasView(first));
  ASSERT_TRUE(grp->hasView(second));
  // Conduct a pointer comparison to make sure that the two views actually
  // refer to the same block of memory
  const double* first_data = grp->getView(first)->getData();
  const double* second_data = grp->getView(second)->getData();
  EXPECT_EQ(first_data, second_data);
}

TEST(sidre_datacollection, create_matset)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  MFEMSidreDataCollection sdc(testName(), &mesh);
  sdc.AssociateMaterialSet("volume_fraction", "matset");

  mfem::H1_FECollection fec(1, mesh.Dimension());
  mfem::FiniteElementSpace fes(&mesh, &fec);

  mfem::GridFunction vol_frac_1(&fes);
  vol_frac_1 = 1.0;

  sdc.RegisterField("volume_fraction_001", &vol_frac_1);

  EXPECT_TRUE(sdc.verifyMeshBlueprint());

  const auto bp_grp = sdc.GetBPGroup();
  EXPECT_TRUE(bp_grp->hasGroup("matsets"));
  EXPECT_TRUE(bp_grp->hasGroup("matsets/matset"));
  checkReferentialEquality(bp_grp,
                           "matsets/matset/volume_fractions/001",
                           "fields/volume_fraction_001/values");
}

TEST(sidre_datacollection, create_matset_multi_fraction)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  MFEMSidreDataCollection sdc(testName(), &mesh);
  sdc.AssociateMaterialSet("volume_fraction", "matset");

  mfem::H1_FECollection fec(1, mesh.Dimension());
  mfem::FiniteElementSpace fes(&mesh, &fec);

  mfem::GridFunction vol_frac_1(&fes);
  vol_frac_1 = 0.7;
  sdc.RegisterField("volume_fraction_001", &vol_frac_1);

  mfem::GridFunction vol_frac_2(&fes);
  vol_frac_2 = 0.3;
  sdc.RegisterField("volume_fraction_002", &vol_frac_2);

  EXPECT_TRUE(sdc.verifyMeshBlueprint());

  const auto bp_grp = sdc.GetBPGroup();
  EXPECT_TRUE(bp_grp->hasGroup("matsets"));
  EXPECT_TRUE(bp_grp->hasGroup("matsets/matset"));
  checkReferentialEquality(bp_grp,
                           "matsets/matset/volume_fractions/001",
                           "fields/volume_fraction_001/values");
  checkReferentialEquality(bp_grp,
                           "matsets/matset/volume_fractions/002",
                           "fields/volume_fraction_002/values");
}

TEST(sidre_datacollection, create_specset)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  MFEMSidreDataCollection sdc(testName(), &mesh);
  mfem::H1_FECollection fec(1, mesh.Dimension());
  mfem::FiniteElementSpace fes(&mesh, &fec);

  // Add a material set so we can associate the specset with it
  sdc.AssociateMaterialSet("volume_fraction", "matset");
  mfem::GridFunction vol_frac_1(&fes);
  vol_frac_1 = 1.0;
  sdc.RegisterField("volume_fraction_001", &vol_frac_1);

  // Then add the species set
  const bool volume_dependent = false;
  sdc.AssociateSpeciesSet("partial_density", "specset", "matset", volume_dependent);

  mfem::GridFunction partial_density_1_1(&fes);
  partial_density_1_1 = 1.0;

  sdc.RegisterField("partial_density_001_001", &partial_density_1_1);

  EXPECT_TRUE(sdc.verifyMeshBlueprint());

  const auto bp_grp = sdc.GetBPGroup();
  EXPECT_TRUE(bp_grp->hasGroup("specsets"));
  EXPECT_TRUE(bp_grp->hasGroup("specsets/specset"));
  EXPECT_TRUE(bp_grp->hasView("specsets/specset/volume_dependent"));
  EXPECT_FALSE(static_cast<axom::int8>(
    bp_grp->getView("specsets/specset/volume_dependent")->getScalar()));
  EXPECT_TRUE(bp_grp->hasView("specsets/specset/matset"));
  EXPECT_EQ(std::string(bp_grp->getView("specsets/specset/matset")->getString()),
            "matset");
  EXPECT_TRUE(bp_grp->hasView("specsets/specset/matset_values/001/001"));
  checkReferentialEquality(bp_grp,
                           "matsets/matset/volume_fractions/001",
                           "fields/volume_fraction_001/values");
  checkReferentialEquality(bp_grp,
                           "specsets/specset/matset_values/001/001",
                           "fields/partial_density_001_001/values");
}

TEST(sidre_datacollection, create_specset_multi_fraction)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  MFEMSidreDataCollection sdc(testName(), &mesh);
  mfem::H1_FECollection fec(1, mesh.Dimension());
  mfem::FiniteElementSpace fes(&mesh, &fec);

  // Add a material set so we can associate the specset with it
  sdc.AssociateMaterialSet("volume_fraction", "matset");
  mfem::GridFunction vol_frac_1(&fes);
  vol_frac_1 = 0.7;
  sdc.RegisterField("volume_fraction_001", &vol_frac_1);
  mfem::GridFunction vol_frac_2(&fes);
  vol_frac_2 = 0.3;
  sdc.RegisterField("volume_fraction_002", &vol_frac_2);

  // Then add the species set
  const bool volume_dependent = false;
  sdc.AssociateSpeciesSet("partial_density", "specset", "matset", volume_dependent);

  mfem::GridFunction partial_density_1_1(&fes);
  partial_density_1_1 = 0.4;
  sdc.RegisterField("partial_density_001_001", &partial_density_1_1);
  mfem::GridFunction partial_density_1_2(&fes);
  partial_density_1_2 = 0.3;
  sdc.RegisterField("partial_density_001_002", &partial_density_1_2);
  mfem::GridFunction partial_density_2_1(&fes);
  partial_density_2_1 = 0.2;
  sdc.RegisterField("partial_density_002_001", &partial_density_2_1);
  mfem::GridFunction partial_density_2_2(&fes);
  partial_density_2_2 = 0.1;
  sdc.RegisterField("partial_density_002_002", &partial_density_2_2);

  EXPECT_TRUE(sdc.verifyMeshBlueprint());

  const auto bp_grp = sdc.GetBPGroup();
  EXPECT_TRUE(bp_grp->hasGroup("specsets"));
  EXPECT_TRUE(bp_grp->hasGroup("specsets/specset"));
  EXPECT_TRUE(bp_grp->hasView("specsets/specset/volume_dependent"));
  EXPECT_FALSE(static_cast<axom::int8>(
    bp_grp->getView("specsets/specset/volume_dependent")->getScalar()));
  EXPECT_TRUE(bp_grp->hasView("specsets/specset/matset"));
  EXPECT_EQ(std::string(bp_grp->getView("specsets/specset/matset")->getString()),
            "matset");
  checkReferentialEquality(bp_grp,
                           "matsets/matset/volume_fractions/001",
                           "fields/volume_fraction_001/values");
  checkReferentialEquality(bp_grp,
                           "matsets/matset/volume_fractions/002",
                           "fields/volume_fraction_002/values");
  checkReferentialEquality(bp_grp,
                           "specsets/specset/matset_values/001/001",
                           "fields/partial_density_001_001/values");
  checkReferentialEquality(bp_grp,
                           "specsets/specset/matset_values/001/002",
                           "fields/partial_density_001_002/values");
  checkReferentialEquality(bp_grp,
                           "specsets/specset/matset_values/002/001",
                           "fields/partial_density_002_001/values");
  checkReferentialEquality(bp_grp,
                           "specsets/specset/matset_values/002/002",
                           "fields/partial_density_002_002/values");
}

TEST(sidre_datacollection, create_material_dependent_field)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  MFEMSidreDataCollection sdc(testName(), &mesh);
  mfem::H1_FECollection fec(1, mesh.Dimension());
  mfem::FiniteElementSpace fes(&mesh, &fec);

  // Add a material set so we can associate the field with it
  sdc.AssociateMaterialSet("volume_fraction", "matset");
  mfem::GridFunction vol_frac_1(&fes);
  vol_frac_1 = 1.0;
  sdc.RegisterField("volume_fraction_001", &vol_frac_1);

  // Then mark the field as material-dependent
  sdc.AssociateMaterialDependentField("density", "matset");

  mfem::GridFunction density_indenpendent(&fes);
  density_indenpendent = 1.0;
  sdc.RegisterField("density", &density_indenpendent);

  mfem::GridFunction density_dependent(&fes);
  density_dependent = 0.2;
  sdc.RegisterField("density_001", &density_dependent);

  EXPECT_TRUE(sdc.verifyMeshBlueprint());

  const auto bp_grp = sdc.GetBPGroup();
  EXPECT_TRUE(bp_grp->hasGroup("fields/density"));
  EXPECT_TRUE(bp_grp->hasGroup("fields/density/matset_values"));
  checkReferentialEquality(bp_grp,
                           "matsets/matset/volume_fractions/001",
                           "fields/volume_fraction_001/values");
  checkReferentialEquality(bp_grp,
                           "fields/density/matset_values/001",
                           "fields/density_001/values");
}

TEST(sidre_datacollection, create_material_dependent_field_multi_fraction)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  MFEMSidreDataCollection sdc(testName(), &mesh);
  mfem::H1_FECollection fec(1, mesh.Dimension());
  mfem::FiniteElementSpace fes(&mesh, &fec);

  // Add a material set so we can associate the field with it
  sdc.AssociateMaterialSet("volume_fraction", "matset");
  mfem::GridFunction vol_frac_1(&fes);
  vol_frac_1 = 0.65;
  sdc.RegisterField("volume_fraction_001", &vol_frac_1);
  mfem::GridFunction vol_frac_2(&fes);
  vol_frac_2 = 0.35;
  sdc.RegisterField("volume_fraction_002", &vol_frac_2);

  // Then mark the field as material-dependent
  sdc.AssociateMaterialDependentField("density", "matset");

  mfem::GridFunction density_indenpendent(&fes);
  density_indenpendent = 1.0;
  sdc.RegisterField("density", &density_indenpendent);

  mfem::GridFunction density_dependent_1(&fes);
  density_dependent_1 = 0.2;
  sdc.RegisterField("density_001", &density_dependent_1);

  mfem::GridFunction density_dependent_2(&fes);
  density_dependent_2 = 0.3;
  sdc.RegisterField("density_002", &density_dependent_2);

  EXPECT_TRUE(sdc.verifyMeshBlueprint());

  const auto bp_grp = sdc.GetBPGroup();
  EXPECT_TRUE(bp_grp->hasGroup("fields/density"));
  EXPECT_TRUE(bp_grp->hasGroup("fields/density/matset_values"));
  checkReferentialEquality(bp_grp,
                           "matsets/matset/volume_fractions/001",
                           "fields/volume_fraction_001/values");
  checkReferentialEquality(bp_grp,
                           "matsets/matset/volume_fractions/002",
                           "fields/volume_fraction_002/values");
  checkReferentialEquality(bp_grp,
                           "fields/density/matset_values/001",
                           "fields/density_001/values");
  checkReferentialEquality(bp_grp,
                           "fields/density/matset_values/002",
                           "fields/density_002/values");
}

  //----------------------------------------------------------------------
  #include "axom/slic/core/SimpleLogger.hpp"
using axom::slic::SimpleLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  SimpleLogger logger;  // create & initialize test logger,

  MPI_Init(&argc, &argv);
  result = RUN_ALL_TESTS();
  MPI_Finalize();

  return result;
}

#endif
