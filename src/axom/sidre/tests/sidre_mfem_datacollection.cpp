// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#ifdef AXOM_USE_MFEM

  #include "mfem.hpp"

  #include "gtest/gtest.h"

  #include "axom/sidre/core/sidre.hpp"
  #include "axom/sidre/core/MFEMSidreDataCollection.hpp"

using axom::sidre::Group;
using axom::sidre::MFEMSidreDataCollection;

const std::string COLL_NAME = "test_collection";
const double EPSILON = 1.0e-6;

TEST(sidre_datacollection, dc_alloc_no_mesh)
{
  MFEMSidreDataCollection sdc(COLL_NAME);
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_alloc_owning_mesh)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  bool owns_mesh = true;
  MFEMSidreDataCollection sdc(COLL_NAME, &mesh, owns_mesh);
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_alloc_nonowning_mesh)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  bool owns_mesh = false;
  MFEMSidreDataCollection sdc(COLL_NAME, &mesh, owns_mesh);
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_register_empty_field)
{
  MFEMSidreDataCollection sdc(COLL_NAME);
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

  MFEMSidreDataCollection sdc(COLL_NAME, &mesh);
  sdc.RegisterField("test_field", &gf);
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_update_state)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  MFEMSidreDataCollection sdc(COLL_NAME, &mesh);

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
  MFEMSidreDataCollection sdc(COLL_NAME, &mesh);

  // The data produced isn't important for this, so just throw it in /tmp
  sdc.SetPrefixPath("/tmp/dc_save_test");
  sdc.Save();

  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_reload)
{
  const std::string field_name = "test_field";
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  mfem::H1_FECollection fec(1, 1);
  mfem::FiniteElementSpace fes(&mesh, &fec);

  // The mesh and field(s) must be owned by Sidre to properly manage data in case of
  // a simulated restart (save -> load)
  bool owns_mesh = true;
  MFEMSidreDataCollection sdc_writer(COLL_NAME, &mesh, owns_mesh);
  mfem::GridFunction gf_write(&fes, nullptr);

  // Register to allocate storage internally, then write to it
  sdc_writer.RegisterField(field_name, &gf_write);

  mfem::ConstantCoefficient three_and_a_half(3.5);
  gf_write.ProjectCoefficient(three_and_a_half);

  EXPECT_TRUE(sdc_writer.verifyMeshBlueprint());

  sdc_writer.SetPrefixPath("/tmp/dc_reload_test");
  sdc_writer.SetCycle(0);
  sdc_writer.PrepareToSave();
  sdc_writer.Save();

  MFEMSidreDataCollection sdc_reader(COLL_NAME, &mesh, owns_mesh);

  sdc_reader.SetPrefixPath("/tmp/dc_reload_test");
  sdc_reader.Load();
  // In this simulated restart a Load deletes the mesh, so it needs to be loaded in again
  sdc_reader.SetMesh(&mesh);

  // Simulate the initial setup process, except this
  // time, the gridfunction data will be overwritten
  mfem::GridFunction gf_read(&fes, nullptr);
  sdc_reader.RegisterField(field_name, &gf_read);

  // Make sure the gridfunction was actually overwritten
  EXPECT_LT(gf_read.ComputeL2Error(three_and_a_half), EPSILON);

  EXPECT_TRUE(sdc_reader.verifyMeshBlueprint());
}

  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
TEST(sidre_datacollection, dc_alloc_owning_parmesh)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  mfem::ParMesh parmesh(MPI_COMM_WORLD, mesh);
  bool owns_mesh = true;
  MFEMSidreDataCollection sdc(COLL_NAME, &parmesh, owns_mesh);
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_alloc_nonowning_parmesh)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  mfem::ParMesh parmesh(MPI_COMM_WORLD, mesh);
  bool owns_mesh = false;
  MFEMSidreDataCollection sdc(COLL_NAME, &parmesh, owns_mesh);
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

  // Replacing data in the mesh with data in the blueprint is not currently
  // supported for parallel meshes, per a comment in createMeshBlueprintAdjacencies

  /*
TEST(sidre_datacollection, dc_par_reload)
{
  const std::string field_name = "test_field";
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  mfem::ParMesh parmesh(MPI_COMM_WORLD, mesh);
  mfem::H1_FECollection fec(1, 1);
  mfem::ParFiniteElementSpace parfes(&parmesh, &fec);

  // The mesh must be owned by Sidre to properly manage data in case of
  // a simulated restart (save -> load)
  bool owns_mesh = true;
  MFEMSidreDataCollection sdc_writer(COLL_NAME, &parmesh, owns_mesh);

  // The mesh and field(s) must be owned by Sidre to properly manage data in case of
  // a simulated restart (save -> load)
  mfem::ParGridFunction gf_write(&parfes, static_cast<double*>(nullptr));

  // Register to allocate storage internally, then write to it
  sdc_writer.RegisterField(field_name, &gf_write);

  mfem::ConstantCoefficient three_and_a_half(3.5);
  gf_write.ProjectCoefficient(three_and_a_half);

  sdc_writer.SetPrefixPath("/tmp/dc_par_reload_test");
  sdc_writer.SetCycle(0);
  sdc_writer.PrepareToSave();
  sdc_writer.Save();

  MFEMSidreDataCollection sdc_reader(COLL_NAME, &parmesh, owns_mesh);
  sdc_reader.SetComm(MPI_COMM_WORLD);
  sdc_reader.SetPrefixPath("/tmp/dc_par_reload_test");
  sdc_reader.Load();
  // In this simulated restart a Load deletes the mesh, so it needs to be loaded in again
  sdc_reader.SetMesh(&parmesh);

  // Simulate the initial setup process, except this
  // time, the gridfunction data will be overwritten
  mfem::ParGridFunction gf_read(&parfes, static_cast<double*>(nullptr));
  sdc_reader.RegisterField(field_name, &gf_read);

  // Make sure the gridfunction was actually overwritten
  EXPECT_LT(gf_read.ComputeL2Error(three_and_a_half), EPSILON);

  EXPECT_TRUE(sdc_reader.verifyMeshBlueprint());
}
*/

    //----------------------------------------------------------------------
    #include "axom/slic/core/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  MPI_Init(&argc, &argv);
  result = RUN_ALL_TESTS();
  MPI_Finalize();

  return result;
}

  #endif

#endif
