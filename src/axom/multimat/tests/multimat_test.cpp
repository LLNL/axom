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
  
  std::vector<int> matcount;

  std::vector<double> volfrac_cellcen_dense;
  std::vector<double> volfrac_matcen_dense;
  std::vector<double> volfrac_cellcen_sparse;
  std::vector<double> volfrac_matcen_sparse;

  std::vector<DataType> cellmat_dense_arr;
  std::vector<DataType> cellmat_sparse_arr;
  std::vector<DataType> matcell_dense_arr;
  std::vector<DataType> matcell_sparse_arr;

  int num_cells, num_mats, stride;
  int nfilled = 0;

  //constructor
  MM_test_data(int n_c = 20, int n_m = 10, int str = 1)
  {
    num_cells = n_c;
    num_mats = n_m;
    stride = str;

    nfilled = 0;
    fillBool_cellcen.resize(num_mats * num_cells, false);
    fillBool_matcen.resize(num_mats * num_cells, false);
    matcount.resize(num_cells, 0);
    for (int c = 0; c < num_cells; ++c) {
      for (int m = 0; m < num_mats; ++m) {
        int dense_cellcen_idx = c * num_mats + m;
        if (dense_cellcen_idx % 3 == 1) {
          int dense_matcen_idx = m * num_cells + c;
          fillBool_cellcen[dense_cellcen_idx] = true;
          fillBool_matcen[dense_matcen_idx] = true;
          nfilled++;
          matcount[c]++;
        }
      }
    }

    //create volfrac array
    int sparse_idx = 0;
    volfrac_cellcen_dense.resize(num_cells*num_mats, 0.0);
    volfrac_matcen_dense.resize(num_cells*num_mats, 0.0);
    volfrac_cellcen_sparse.resize(nfilled, 0.0);
    volfrac_matcen_sparse.resize(nfilled, 0.0);
    for (auto i = 0; i < num_cells; ++i)
    {
      for (auto m = 0; m < num_mats; ++m) {
        if (fillBool_cellcen[i*num_mats + m]) {
          volfrac_cellcen_dense[i*num_mats + m] = 1.0 / (double)matcount[i];
          volfrac_matcen_dense[m*num_cells + i] = 1.0 / (double)matcount[i];
          volfrac_cellcen_sparse[sparse_idx++] = 1.0 / (double)matcount[i];
        }
      }
    }
    sparse_idx = 0;
    for (int m = 0; m < num_mats; ++m) {
      for (int c = 0; c < num_cells; ++c) {
        int dense_cellcen_idx = c * num_mats + m;
        if (fillBool_cellcen[dense_cellcen_idx]) {
          volfrac_matcen_sparse[sparse_idx++] = 1.0 / (double)matcount[c];
        }
      }
    }

    //fill in the data
    cellmat_dense_arr.resize(num_mats * num_cells * stride);
    cellmat_sparse_arr.resize(nfilled * stride);
    matcell_dense_arr.resize(num_mats * num_cells * stride);
    matcell_sparse_arr.resize(nfilled * stride);
    sparse_idx = 0;
    for (int c = 0; c < num_cells; ++c) {
      for (int m = 0; m < num_mats; ++m) {
        int dense_cellcen_idx = c * num_mats + m;
        int dense_matcen_idx = m * num_cells + c;

        if (fillBool_cellcen[dense_cellcen_idx]) {
          for (int s = 0; s < stride; ++s) {
            DataType val = static_cast<DataType>(c * 1000.0 + m * 1.0 + s * 0.01);
            cellmat_sparse_arr[sparse_idx*stride + s] = val;
            cellmat_dense_arr[dense_cellcen_idx*stride + s] = val;
          }
          ++sparse_idx;
        }
        if (fillBool_cellcen[dense_cellcen_idx]) {
          for (int s = 0; s < stride; ++s) {
            DataType val = static_cast<DataType>(c * 1000.0 + m * 1.0 + s * 0.01);
            matcell_dense_arr[dense_matcen_idx*stride + s] = val;
          }
        }
      }
    }
    sparse_idx = 0;
    for (int m = 0; m < num_mats; ++m) {
      for (int c = 0; c < num_cells; ++c) {
        int dense_cellcen_idx = c * num_mats + m;
        if (fillBool_cellcen[dense_cellcen_idx]) {
          for (int s = 0; s < stride; ++s) {
            DataType val = cellmat_dense_arr[dense_cellcen_idx*stride + s];
            matcell_sparse_arr[sparse_idx*stride + s] = val;
          }
          ++sparse_idx;
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
  EXPECT_EQ(data.stride, map.stride());
  
  //check via findValue(...)
  int sparse_idx = 0;
  for (int ci = 0; ci < data.num_cells; ++ci)
  {
    for (int mi = 0; mi < data.num_mats; ++mi)
    {
      for (int s = 0; s < data.stride; ++s)
      {
        double *d;
        if (mm.getDataLayout() == DataLayout::CELL_CENTRIC)
          d = map.findValue(ci, mi, s);
        else
          d = map.findValue(mi, ci, s);
        int dense_idx = ci * data.num_mats + mi;

        if (mm.getSparcityLayout() == SparcityLayout::DENSE) {
          EXPECT_EQ(*d, data.cellmat_dense_arr[dense_idx*data.stride + s]);
        }
        else {
          if (data.fillBool_cellcen[dense_idx]) {
            EXPECT_NE(d, nullptr);
            EXPECT_EQ(*d, data.cellmat_dense_arr[dense_idx*data.stride + s]);
            if(s == data.stride-1)
              sparse_idx += 1;
          }
          else {
            EXPECT_EQ(d, nullptr);
          }
        }
      }
    }
  }
}

template<typename T>
MultiMat* newMM(MM_test_data<T>& data, DataLayout layout_used, SparcityLayout sparcity_used, std::string& array_name)
{
  MultiMat* mm_ptr = new MultiMat(layout_used, sparcity_used);
  MultiMat& mm = *mm_ptr;
  mm.setNumberOfCell(data.num_cells);
  mm.setNumberOfMat(data.num_mats);

  if (layout_used == DataLayout::CELL_CENTRIC)
    mm.setCellMatRel(data.fillBool_cellcen);
  else
    mm.setCellMatRel(data.fillBool_matcen);

  if (layout_used == DataLayout::CELL_CENTRIC) {
    if (sparcity_used == SparcityLayout::DENSE)
      mm.setVolfracField(data.volfrac_cellcen_dense.data());
    else
      mm.setVolfracField(data.volfrac_cellcen_sparse.data());
  }
  else {
    if (sparcity_used == SparcityLayout::DENSE)
      mm.setVolfracField(data.volfrac_matcen_dense.data());
    else
      mm.setVolfracField(data.volfrac_matcen_sparse.data());
  }

  EXPECT_TRUE(mm.isValid());

  if (layout_used == DataLayout::CELL_CENTRIC) {
    if (sparcity_used == SparcityLayout::DENSE)
      mm.addField(array_name, FieldMapping::PER_CELL_MAT, data.cellmat_dense_arr.data(), data.stride);
    else
      mm.addField(array_name, FieldMapping::PER_CELL_MAT, data.cellmat_sparse_arr.data(), data.stride);
  }
  else {
    if (sparcity_used == SparcityLayout::DENSE)
      mm.addField(array_name, FieldMapping::PER_CELL_MAT, data.matcell_dense_arr.data(), data.stride);
    else
      mm.addField(array_name, FieldMapping::PER_CELL_MAT, data.matcell_sparse_arr.data(), data.stride);
  }

  return mm_ptr;
}


TEST(multimat, construct_multimat_1_array)
{
  const int num_cells = 20;
  const int num_mats = 10;
  const int stride_val = 4;
  MM_test_data<double> data(num_cells, num_mats, stride_val);

  std::vector<DataLayout> data_layouts = { DataLayout::CELL_CENTRIC, DataLayout::MAT_CENTRIC };
  std::vector<SparcityLayout> sparcity_layouts = { SparcityLayout::DENSE, SparcityLayout::SPARSE };

  std::string array_name = "Array 1";

  for (auto layout_used : data_layouts) {
    for (auto sparcity_used : sparcity_layouts) {
      SLIC_INFO("--------------------\nConstructing MultiMat object...");
      MultiMat* mm_ptr = newMM(data, layout_used, sparcity_used, array_name);
      MultiMat& mm = *mm_ptr;

      EXPECT_TRUE(mm.isValid(true));
      EXPECT_EQ(mm.getDataLayout(), layout_used);
      EXPECT_EQ(mm.getSparcityLayout(), sparcity_used);
      EXPECT_EQ(mm.getNumberOfCells(), data.num_cells);
      EXPECT_EQ(mm.getNumberOfMaterials(), data.num_mats);

      check_values<double>(mm, array_name, data);

      for (auto layout_to_convert : data_layouts)
      {
        for (auto sparcity_to_convert : sparcity_layouts)
        {
          SLIC_INFO("Making a copy...");
          //Make Copies to test conversion
          MultiMat mm_c(mm);
          EXPECT_TRUE(mm_c.isValid(true));
          EXPECT_EQ(mm_c.getDataLayout(), layout_used);
          EXPECT_EQ(mm_c.getSparcityLayout(), sparcity_used);
          EXPECT_EQ(mm_c.getNumberOfCells(), data.num_cells);
          EXPECT_EQ(mm_c.getNumberOfMaterials(), data.num_mats);

          check_values<double>(mm_c, array_name, data);

          SLIC_INFO("Converting layout...");
          mm_c.convertLayout(layout_to_convert, sparcity_to_convert);
          SLIC_INFO("Layout converted from " << mm.getDataLayoutAsString()
            << " and " << mm.getSparcityLayoutAsString() << 
            "\n                          to " << mm_c.getDataLayoutAsString() 
            << " and " << mm_c.getSparcityLayoutAsString() <<"...");

          EXPECT_TRUE(mm_c.isValid(true));
          EXPECT_EQ(mm_c.getDataLayout(), layout_to_convert);
          EXPECT_EQ(mm_c.getSparcityLayout(), sparcity_to_convert);
          EXPECT_EQ(mm_c.getNumberOfCells(), data.num_cells);
          EXPECT_EQ(mm_c.getNumberOfMaterials(), data.num_mats);
          check_values<double>(mm_c, array_name, data);

        }
      }

      delete mm_ptr;
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
