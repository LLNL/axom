// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#ifdef AXOM_USE_MFEM

#include "mfem.hpp"

#include "gtest/gtest.h"

#include "axom/sidre/core/sidre.hpp"
#include "axom/sidre/datacollection/SidreDataCollection.hpp"

using axom::sidre::Group;

const std::string COLL_NAME = "test_collection";

TEST(sidre_datacollection, dc_alloc_no_mesh)
{
    mfem::SidreDataCollection sdc(COLL_NAME);
    EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_alloc_owning_mesh)
{
    // 1D mesh divided into 10 segments
    mfem::Mesh mesh(10);
    bool owns_mesh = true;
    mfem::SidreDataCollection sdc(COLL_NAME, &mesh, owns_mesh);
    EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_alloc_nonowning_mesh)
{
    // 1D mesh divided into 10 segments
    mfem::Mesh mesh(10);
    bool owns_mesh = false;
    mfem::SidreDataCollection sdc(COLL_NAME, &mesh, owns_mesh);
    EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_register_empty_field)
{
    mfem::SidreDataCollection sdc(COLL_NAME);
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
    
    mfem::SidreDataCollection sdc(COLL_NAME, &mesh);
    sdc.RegisterField("test_field", &gf);
    EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_update_state)
{
    // 1D mesh divided into 10 segments
    mfem::Mesh mesh(10);
    mfem::SidreDataCollection sdc(COLL_NAME, &mesh);

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
    mfem::SidreDataCollection sdc(COLL_NAME, &mesh);
    
    // The data produced isn't important for this, so just throw it in /tmp
    sdc.SetPrefixPath("/tmp/dc_save_test");
    sdc.Save();

    EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_reload)
{
    // 1D mesh divided into 10 segments
    mfem::Mesh mesh(10);
    // The mesh must be owned by Sidre to properly manage data in case of
    // a simulated restart (save -> load)
    bool owns_mesh = true;
    mfem::SidreDataCollection sdc_writer(COLL_NAME, &mesh, owns_mesh);

    sdc_writer.SetPrefixPath("/tmp/dc_reload_test");
    sdc_writer.SetCycle(0);
    sdc_writer.PrepareToSave();
    sdc_writer.Save();
    
    mfem::SidreDataCollection sdc_reader(COLL_NAME);
    sdc_reader.SetPrefixPath("/tmp/dc_reload_test");
    sdc_reader.Load();

    EXPECT_TRUE(sdc_reader.verifyMeshBlueprint());
}



#endif
