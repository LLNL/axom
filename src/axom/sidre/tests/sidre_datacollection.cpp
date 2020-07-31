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

TEST(sidre_datacollection, dc_alloc_no_mesh)
{
    mfem::SidreDataCollection sdc("test_collection");

    Group* bp_group = sdc.GetBPGroup();
    Group* bp_idx_group = sdc.GetBPIndexGroup();

    // Should be empty at initialization when no mesh passed
    EXPECT_EQ(bp_group->getNumGroups(), 0);
    EXPECT_EQ(bp_idx_group->getNumGroups(), 0);
}

TEST(sidre_datacollection, dc_alloc_mesh)
{
    mfem::Mesh mesh;
    mfem::SidreDataCollection sdc("test_collection", &mesh);

    Group* bp_group = sdc.GetBPGroup();
    Group* bp_idx_group = sdc.GetBPIndexGroup();

    // Should be empty at initialization when no mesh passed
    EXPECT_EQ(bp_group->getNumGroups(), 4);
    EXPECT_EQ(bp_idx_group->getNumGroups(), 4);
}

#endif
