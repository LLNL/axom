/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

#include "gtest/gtest.h"

#include <vector>

#include "sidre/sidre.hpp"

using asctoolkit::sidre::DataBuffer;
using asctoolkit::sidre::DataGroup;
using asctoolkit::sidre::DataStore;
using asctoolkit::sidre::DataView;

using namespace conduit;

//------------------------------------------------------------------------------

TEST(sidre_external, simple_arrays)
{
    DataStore* ds   = new DataStore();
    DataGroup* root = ds->getRoot();

    const int len = 11;

    int*    idata = new int[len]; 
    double* ddata = new double[len];

    for (int ii = 0; ii < len; ++ii) {
       idata[ii] = ii;
       ddata[ii] = idata[ii] * 2.0;
    }

    (void) root->createExternalView("idata", idata, 
                                     DataType::c_int(len));
    (void) root->createExternalView("ddata", ddata, 
                                     DataType::c_double(len));
    EXPECT_EQ(root->getNumViews(), 2u);

    root->getView("idata")->getNode().print_detailed();  
    root->getView("ddata")->getNode().print_detailed();  

    int* idata_chk = root->getView("idata")->getNode().as_int_ptr();
    for (int ii = 0; ii < len; ++ii) {
       EXPECT_EQ(idata_chk[ii], idata[ii]);
    }

    double* ddata_chk = root->getView("ddata")->getNode().as_double_ptr();
    for (int ii = 0; ii < len; ++ii) {
       EXPECT_EQ(ddata_chk[ii], ddata[ii]);
    }

    delete ds;
    delete [] idata;
    delete [] ddata;
}

//------------------------------------------------------------------------------

#if 0  // RDH -- Still think through how to handle this...
TEST(sidre_external, save_restore_simple_arrays)
{
    DataStore* ds   = new DataStore();
    DataGroup* root = ds->getRoot();

    const int len = 11;

    int*    idata = new int[len];
    double* ddata = new double[len];

    for (int ii = 0; ii < len; ++ii) {
       idata[ii] = ii;
       ddata[ii] = idata[ii] * 2.0;
    }

    (void) root->createExternalView("idata", idata,
                                     DataType::c_int(len));
    (void) root->createExternalView("ddata", ddata,
                                     DataType::c_double(len));

#if 0
    root->getView("idata")->getNode().print_detailed();
    root->getView("ddata")->getNode().print_detailed();
#endif

    ds->getRoot()->save("out_sidre_external_save_restore_simple_arrays", "conduit");
    ds->print();


    DataStore *ds2 = new DataStore();

    ds2->getRoot()->load("out_sidre_external_save_restore_simple_arrays","conduit");

    ds2->print();

    DataGroup* root2 = ds2->getRoot();

    int* idata_chk = root2->getView("idata")->getNode().as_int_ptr();
    for (int ii = 0; ii < len; ++ii) {
       EXPECT_EQ(idata_chk[ii], idata[ii]);
    }

    double* ddata_chk = root2->getView("ddata")->getNode().as_double_ptr();
    for (int ii = 0; ii < len; ++ii) {
       EXPECT_EQ(ddata_chk[ii], ddata[ii]);
    }

    delete ds;
    delete [] idata;
    delete [] ddata;
}
#endif
