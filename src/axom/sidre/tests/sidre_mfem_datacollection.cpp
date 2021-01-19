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

  mfem::ConstantCoefficient one(1.0);
  gf_write.ProjectCoefficient(one);

  EXPECT_TRUE(sdc_writer.verifyMeshBlueprint());

  sdc_writer.SetPrefixPath("/tmp/dc_reload_test");
  sdc_writer.SetCycle(0);
  sdc_writer.PrepareToSave();
  sdc_writer.Save();

  MFEMSidreDataCollection sdc_reader(COLL_NAME);
  sdc_reader.SetPrefixPath("/tmp/dc_reload_test");
  sdc_reader.Load();

  EXPECT_TRUE(sdc_reader.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_qf_reload)
{
  //Set up a small mesh and a couple of grid function on that mesh
  mfem::Mesh *mesh = new mfem::Mesh(2, 3, mfem::Element::QUADRILATERAL, 0, 2.0, 3.0);
  mfem::FiniteElementCollection *fec = new mfem::LinearFECollection;
  mfem::FiniteElementSpace *fespace = new mfem::FiniteElementSpace(mesh, fec);

  const int intOrder = 3;
  const int qs_vdim = 1;
  const int qv_vdim = 2;

  mfem::QuadratureSpace *qspace = new mfem::QuadratureSpace(mesh, intOrder);
  // We want Sidre to allocate the data for us and to own the internal data.
  // If we don't do the below then the data collection doesn't save off the data
  // structure.
  mfem::QuadratureFunction *qs = new mfem::QuadratureFunction(qspace, nullptr, qs_vdim);
  mfem::QuadratureFunction *qv = new mfem::QuadratureFunction(qspace, nullptr, qv_vdim);

  int Nq = qs->Size();

  // The mesh and field(s) must be owned by Sidre to properly manage data in case of
  // a simulated restart (save -> load)
  bool owns_mesh = true;
  MFEMSidreDataCollection sdc_writer(COLL_NAME, mesh, owns_mesh);

  // sdc_writer owns the quadrature function fields and the underlying data
  sdc_writer.RegisterQField("qs", qs);
  sdc_writer.RegisterQField("qv", qv);

  // The data needs to be instantiated before we save it off
  for (int i = 0; i < Nq; ++i)
  {
      (*qs)(i) = double(i);
      (*qv)(2*i+0) = double(i);
      (*qv)(2*i+1) = double(Nq - i - 1);
  }

  sdc_writer.SetPrefixPath("/tmp/dc_qf_test");
  sdc_writer.SetCycle(5);
  sdc_writer.SetTime(8.0);
  sdc_writer.Save();

  MFEMSidreDataCollection sdc_reader(COLL_NAME);
  sdc_reader.SetPrefixPath("/tmp/dc_qf_test");
  sdc_reader.Load(sdc_writer.GetCycle());

  const int order_qs = sdc_reader.GetQFieldOrder("qs");
  const int vdim_qs = sdc_reader.GetQFieldVDim("qs");
  double* data_qs = sdc_reader.GetQFieldData("qs");

  const int order_qv = sdc_reader.GetQFieldOrder("qv");
  const int vdim_qv = sdc_reader.GetQFieldVDim("qv");

  // order_qs should also equal order_qv in this trivial case
  EXPECT_TRUE(order_qs == intOrder);
  EXPECT_TRUE(order_qv == intOrder);

  EXPECT_TRUE(vdim_qs == qs_vdim);
  EXPECT_TRUE(vdim_qv == qv_vdim);

  mfem::QuadratureSpace qspace_new(mesh, order_qs);

  // Just showing different ways in which you could restore a QuadratureFunction
  mfem::QuadratureFunction qs_new(&qspace_new, data_qs, vdim_qs);
  mfem::QuadratureFunction qv_new(&qspace_new, sdc_reader.GetQFieldData("qv"), sdc_reader.GetQFieldVDim("qv"));

  qs_new -= *qs;
  qv_new -= *qv;

  EXPECT_TRUE(qs_new.Norml2() < 1e-15);
  EXPECT_TRUE(qv_new.Norml2() < 1e-15);

  delete fespace;
  delete qs; delete qv; delete qspace;
  delete mesh;

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

TEST(sidre_datacollection, dc_par_reload)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  mfem::ParMesh parmesh(MPI_COMM_WORLD, mesh);
  // The mesh must be owned by Sidre to properly manage data in case of
  // a simulated restart (save -> load)
  bool owns_mesh = true;
  MFEMSidreDataCollection sdc_writer(COLL_NAME, &parmesh, owns_mesh);

  sdc_writer.SetPrefixPath("/tmp/dc_par_reload_test");
  sdc_writer.SetCycle(0);
  sdc_writer.PrepareToSave();
  sdc_writer.Save();

  MFEMSidreDataCollection sdc_reader(COLL_NAME);
  sdc_reader.SetComm(MPI_COMM_WORLD);
  sdc_reader.SetPrefixPath("/tmp/dc_par_reload_test");
  sdc_reader.Load();

  EXPECT_TRUE(sdc_reader.verifyMeshBlueprint());
}

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
