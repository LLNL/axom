/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


/**
 * \file multimat_test.cpp
 *
 * \brief Unit tests for MultiMat class
 */

#include <iterator>
#include "gtest/gtest.h"

#include "multimat/multimat.hpp"

#include "slic/slic.hpp"
#include "slic/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

using namespace axom::multimat;

TEST(multimat,construct_empty_multimat_obj)
{
  MultiMat mm;

  EXPECT_TRUE(mm.isValid(true));
}

template<typename DataType>
struct MM_test_data {
  std::vector<bool> fillBool_cellcen;
  std::vector<bool> fillBool_matcen;

  std::vector<DataType> cellmat_dense_arr;
  std::vector<DataType> cellmat_sparse_arr;
  std::vector<DataType> matcell_dense_arr;
  std::vector<DataType> matcell_sparse_arr;

  int num_cells, num_mats;
  int nfilled = 0;

  //constructor
  MM_test_data(int n_c = 20, int n_m = 10)
  {
    num_cells = n_c;
    num_mats = n_m;
    
    nfilled = 0;
    fillBool_cellcen.resize(num_mats * num_cells, false);
    fillBool_matcen.resize(num_mats * num_cells, false);
    for (int c = 0; c < num_cells; ++c) {
      for (int m = 0; m < num_mats; ++m) {
        int dense_cellcen_idx = c * num_mats + m;
        if (dense_cellcen_idx % 3 == 1) {
          int dense_matcen_idx = m * num_cells + c;
          fillBool_cellcen[dense_cellcen_idx] = true;
          fillBool_matcen[dense_matcen_idx] = true;
          nfilled++;
        }
      }
    }

    //fill in the data
    cellmat_dense_arr.resize(num_mats * num_cells);
    cellmat_sparse_arr.resize(nfilled);
    matcell_dense_arr.resize(num_mats * num_cells);
    matcell_sparse_arr.resize(nfilled);
    int sparse_idx = 0;
    for (int c = 0; c < num_cells; ++c) {
      for (int m = 0; m < num_mats; ++m) {
        int dense_cellcen_idx = c * num_mats + m;
        int dense_matcen_idx = m * num_cells + c;
        DataType val = static_cast<DataType>(c * 1000.0 + m * 1.0);
        if (fillBool_cellcen[dense_cellcen_idx]) {
          cellmat_sparse_arr[sparse_idx++] = val;
          cellmat_dense_arr[dense_cellcen_idx] = val;
        }
        if (fillBool_cellcen[dense_cellcen_idx]) {
          matcell_dense_arr[dense_matcen_idx] = val;
        }
      }
    }
    sparse_idx = 0;
    for (int m = 0; m < num_mats; ++m) {
      for (int c = 0; c < num_cells; ++c) {
        int dense_cellcen_idx = c * num_mats + m;
        DataType val = cellmat_dense_arr[dense_cellcen_idx];
        if (fillBool_cellcen[dense_cellcen_idx]) {
          matcell_sparse_arr[sparse_idx++] = val;
        }
      }
    }

  }
};


template<typename DataType>
void check_values(MultiMat& mm, std::string arr_name, MM_test_data<DataType>& data)
{
  MultiMat::Field2D<DataType>& map = mm.get2dField<DataType>(arr_name);
  EXPECT_TRUE(map.isValid());

  //check via findValue(...)
  int sparse_idx = 0;
  for (int ci = 0; ci < data.num_cells; ++ci)
  {
    for (int mi = 0; mi < data.num_mats; ++mi)
    {
      double *d;
      if(mm.getDataLayout()==DataLayout::CELL_CENTRIC) 
        d = map.findValue(ci, mi);
      else
        d = map.findValue(mi, ci);
      int dense_idx = ci * data.num_mats + mi;

      if (mm.getSparcityLayout() == SparcityLayout::DENSE) {
        EXPECT_EQ(*d, data.cellmat_dense_arr[dense_idx]);
      }
      else {
        if (data.fillBool_cellcen[dense_idx]) {
          EXPECT_NE(d, nullptr);
          if (d) EXPECT_EQ(*d, data.cellmat_dense_arr[dense_idx]);
          sparse_idx += 1;
        }
        else {
          EXPECT_EQ(d, nullptr);
        }
      }
    }
  }
}


TEST(multimat, construct_multimat_1_array)
{
  const int num_cells = 20;
  const int num_mats = 10;
  MM_test_data<double> data(20, 10);

  std::vector<DataLayout> data_layouts = { DataLayout::CELL_CENTRIC, DataLayout::MAT_CENTRIC };
  std::vector<SparcityLayout> sparcity_layouts = { SparcityLayout::DENSE, SparcityLayout::SPARSE };

  for (auto layout_used : data_layouts) {
    for (auto sparcity_used : sparcity_layouts) {
      SLIC_INFO("Constructing MultiMat object...");
      MultiMat mm(layout_used, sparcity_used);
      mm.setNumberOfCell(data.num_cells);
      mm.setNumberOfMat(data.num_mats);
      if (layout_used == DataLayout::CELL_CENTRIC) 
        mm.setCellMatRel(data.fillBool_cellcen);
      else
        mm.setCellMatRel(data.fillBool_matcen);

      EXPECT_TRUE(mm.isValid());

      std::string array_name = "Array 1";

      if (layout_used == DataLayout::CELL_CENTRIC) {
        if(sparcity_used == SparcityLayout::DENSE)
          mm.newFieldArray(array_name, FieldMapping::PER_CELL_MAT, data.cellmat_dense_arr.data());
        else
          mm.newFieldArray(array_name, FieldMapping::PER_CELL_MAT, data.cellmat_sparse_arr.data());
      }
      else {
        if (sparcity_used == SparcityLayout::DENSE)
          mm.newFieldArray(array_name, FieldMapping::PER_CELL_MAT, data.matcell_dense_arr.data());
        else
          mm.newFieldArray(array_name, FieldMapping::PER_CELL_MAT, data.matcell_sparse_arr.data());
      }

      EXPECT_TRUE(mm.isValid(true));
      EXPECT_EQ(mm.getDataLayout(), layout_used);
      EXPECT_EQ(mm.getSparcityLayout(), sparcity_used);
      EXPECT_EQ(mm.getNumberOfCells(), data.num_cells);
      EXPECT_EQ(mm.getNumberOfMaterials(), data.num_mats);

      check_values<double>(mm, array_name, data);
    }
  }

}


//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
#ifdef AXOM_DEBUG
  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
#endif

  UnitTestLogger logger;  // create & initialize test logger,
  axom::slic::setLoggingMsgLevel( axom::slic::message::Info );

  int result = RUN_ALL_TESTS();

  return result;
}
