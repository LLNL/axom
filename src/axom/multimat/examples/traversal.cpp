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
 * \file traversal.cpp
 *
 * \brief Example of traversing data using various methods in MultiMat
 */

#include "axom/multimat/multimat.hpp"

#include <ctime>

using namespace std;
using namespace axom::multimat;


std::clock_t start_time;
void start_timer() {
  start_time = std::clock();
}
double end_timer() {
  double duration = (std::clock() - start_time) / (double)CLOCKS_PER_SEC;
  return duration;
}


void various_traversal_methods(int nmats, int ncells, int ncomp, bool use_sparse) {

  MultiMat mm(DataLayout::CELL_CENTRIC, use_sparse ? SparcityLayout::SPARSE : SparcityLayout::DENSE);
  
  int nfilled = 0;
  std::vector<bool> cellMatRel(nmats * ncells, false);
  for (int i = 0; i < nmats*ncells; i++) {
    if (i % 3 == 1) {
      cellMatRel[i] = true;
      nfilled++;
    }
  }

  mm.setNumberOfMat(nmats);
  mm.setNumberOfCell(ncells);
  mm.setCellMatRel(cellMatRel);

  //create the std::vector data for the field arrays
  std::vector<double> cell_arr1(ncells*ncomp);
  double c_sum = 0;
  for (int i = 0; i < ncells; ++i) {
    for (int comp = 0; comp < ncomp; ++comp) {
      cell_arr1[i*ncomp + comp] = (double)i * 2.0 + comp * 0.01;
      c_sum += cell_arr1[i*ncomp + comp];
    }
  }
    
  std::vector<double> cellmat_arr;
  cellmat_arr.resize( (use_sparse ? nfilled : nmats * ncells) *ncomp);
  double x_sum = 0;
  for (unsigned int i = 0; i < cellmat_arr.size()/ncomp; i++) {
    if (use_sparse || cellMatRel[i]) {
      for (int comp = 0; comp < ncomp; ++comp) {
        cellmat_arr[i*ncomp + comp] = (double)i * 1.1 + comp * 0.001;
        x_sum += cellmat_arr[i*ncomp + comp];
      }
    }
  }

  //create volfrac array
  std::vector<double> volfrac_arr(ncells*nmats, 0);
  for (auto i = 0; i < ncells; ++i)
  {
    int matcount = 0;
    for (auto m = 0; m < nmats; ++m) {
      if (cellMatRel[i*nmats + m])
        matcount+=1;
    }
    for (auto m = 0; m < nmats; ++m) {
      if (cellMatRel[i*nmats + m])
        volfrac_arr[i*nmats + m] = 1.0 / (double)matcount;
    }
  }

  mm.setVolfracField(volfrac_arr.data());
  mm.addField("Cell Array"   , FieldMapping::PER_CELL    , &cell_arr1[0], ncomp);
  mm.addField("CellMat Array", FieldMapping::PER_CELL_MAT, &cellmat_arr[0], ncomp);

  //convert layout
  mm.convertLayoutToSparse();

  double sum = 0;

  //Different accessing methods ...
  
  //get the volfrac field
  MultiMat::Field2D<double>& volfrac_map = mm.getVolfracField();
  MultiMat::Field2D<double>& volfrac_map2 = mm.get2dField<double>("Volfrac");
  assert(&volfrac_map == &volfrac_map2);
  
  // --------- returning SLAM map and submap -----------
  printf("\nAccess from SLAM map (and submap)\n");
  sum = 0;
  start_timer();
  {
    MultiMat::Field1D<double>& map = mm.get1dField<double>("Cell Array");
    assert(ncomp == map.stride());
    assert(ncomp == map.numComp());
    for (int i = 0; i < mm.getNumberOfCells(); i++) {
      for (int comp = 0; comp < map.numComp(); comp++)
      {
        double val = map(i, comp);                               //<----
        assert(val == map[ i*map.numComp() + comp] );    //1d bracket access
        sum += val;
      }
    }
  }
  cout << end_timer() << "\t";
  assert(c_sum == sum);

  sum = 0;
  start_timer();
  {
    MultiMat::Field2D<double>& map2d = mm.get2dField<double>("CellMat Array");
    assert(ncomp == map2d.stride());
    assert(ncomp == map2d.numComp());
    for (int i = 0; i < map2d.firstSetSize(); i++)
    {
      MultiMat::IdSet rel_set = map2d.indexSet(i);
      MultiMat::SubField<double> submap = map2d(i);
      assert(rel_set.size() == submap.size());
      for (int k = 0; k < submap.size(); k++)
      {
        int idx = submap.index(k);  //mat id
        int idx2 = rel_set[k]; //another way to get mat id
        assert(idx == idx2); 

        for (int c = 0; c < submap.numComp(); ++c) {
          double val = submap.value(k, c);                 //<----------
          assert( val == submap(k,c) );                    //operator () access
          assert( val == submap[ k*submap.numComp()+c ] );  //1d bracket access
          sum += val;         
        }
      }
    }
  }
  cout << end_timer() << "\t";
  assert(x_sum == sum);
  

  // ------- Dense Access ----------
  printf("\nDense Access via map\n-\t");
  sum = 0;
  start_timer();
  {
    MultiMat::Field2D<double>& map = mm.get2dField<double>("CellMat Array");
    
    for (int i = 0; i < mm.getNumberOfCells(); i++) {
      for (int m = 0; m < mm.getNumberOfMaterials(); m++) {
        for (int c = 0; c < map.numComp(); ++c) 
        {
          double* valptr = map.findValue(i, m, c); //<---- contains a hidden for-loop for sparse layouts
          if (valptr) 
            sum += *valptr;           
        }
      }
    }
  }
  cout << end_timer() << "\t";
  assert(x_sum == sum);
   

  if (mm.getSparcityLayout() == SparcityLayout::SPARSE)
  {
    // ------------ return index set --------------
    printf("\nAccess by Map with indexing set\n-\t");
    sum = 0;
    start_timer();
    {
      //Note you can achieve the same thing by using a Map pointer to point to a bivariateMap object
      //MultiMat::Field1D<double>& map = mm.get1dField<double>("CellMat Array");

      MultiMat::Field2D<double>& map = mm.get2dField<double>("CellMat Array");

      for (int i = 0; i < mm.getNumberOfCells(); i++)
      {
        MultiMat::IdSet setOfMaterialsInThisCell = mm.getMatInCell(i); //the materials (by id) in this cell
        MultiMat::IndexSet indexSet = mm.getIndexingSetOfCell(i); //the indices into the maps
        assert(setOfMaterialsInThisCell.size() == indexSet.size());
        for (int j = 0; j < indexSet.size(); j++)
        {
          int mat_id = setOfMaterialsInThisCell.at(j);
          for (int comp = 0; comp < map.numComp(); ++comp) {
            double val = map[indexSet[j]*map.numComp() + comp];   //<-----
            //assert(val == map(indexSet[j], comp)); //if 1dMap is used, this is also valid
            sum += val;
          }
        }
      }
    }
    cout << end_timer() << "\t";
    assert(x_sum == sum);
  }


  // ---------- using iterator with Map and Submap -------------
  printf("\nWith Map (and Submap) iterators\n");
  sum = 0;
  start_timer();
  {
    MultiMat::Field1D<double>& map = mm.get1dField<double>("Cell Array");
    for (MultiMat::Field1D<double>::iterator iter = map.begin(); iter != map.end(); iter++)
    {
      for(int comp=0; comp<iter.numComp(); ++comp)
        sum += iter(comp);              //<----
    }
  }
  cout << end_timer() << "\t";
  assert(c_sum == sum);


  sum = 0;
  start_timer();
  {
    MultiMat::Field2D<double>& map2d = mm.get2dField<double>("CellMat Array");
    for (int i = 0; i < map2d.firstSetSize(); i++) 
    {
      MultiMat::SubField<double> submap = map2d(i);
      for (auto iter = submap.begin(); iter != submap.end(); iter++) 
      {
        for (int comp = 0; comp < map2d.numComp(); ++comp)
        {
          sum += iter(comp);                  //<----
          assert(iter(comp) == iter.value(comp));  //another way to get the value
          int idx = iter.index();
        }
        //assert(*iter, iter.value());   //2 ways to get the first component value
      }
    }
  }
  cout << end_timer() << "\t";
  assert(x_sum == sum);


  // ---------- iterator for BivariateMap ------------
  printf("\nWith Map iterators - begin(i) and end(i)\n-\t");
  sum = 0;
  start_timer();
  {
    MultiMat::Field2D<double>& map2d = mm.get2dField<double>("CellMat Array");
    for (int i = 0; i < mm.getNumberOfCells();/* map2d.firstSetSize(); */i++)
    {
      for (MultiMat::SubField<double>::SubMapIterator iter = map2d.begin(i); iter != map2d.end(i); iter++)
      {
        for (int comp = 0; comp < map2d.numComp(); ++comp)
        {
          sum += iter(comp);          //<----
          assert(iter(comp) == iter.value(comp));  //another way to get the value
          int idx = iter.index();
        }
        assert(iter(0) == *iter); //2 ways to get the first component value
        assert(iter(0) == iter.value()); 
      }
    }
  }
  cout << end_timer() << "\t";
  assert(x_sum == sum);



  // ---------- range-based for-loop, only works if there is 1 component ------------
  if( ncomp == 1 )
  {
    printf("\nWith range-based for-loop of Map \n-\t");
    sum = 0;
    start_timer();
    
    MultiMat::Field2D<double>& map2d = mm.get2dField<double>("CellMat Array");
    for (int i = 0; i < mm.getNumberOfCells(); i++)
    {
      MultiMat::SubField<double> submap = map2d(i);
      for (double val : submap)
      {
        sum += val;          //<----
      }
    }
    
    cout << end_timer() << "\t";
    assert(x_sum == sum);
  }


  // ---------- flat iterator ------------
  printf("\nWith BivariateMap flat iterators \n-\t");
  sum = 0;
  start_timer();
  {
    MultiMat::Field2D<double>& map2d = mm.get2dField<double>("CellMat Array");
    for (auto iter = map2d.begin(); iter != map2d.end(); ++iter)
    {
      int cell_id = iter.firstIndex();
      int mat_id = iter.secondIndex();
      for (int comp = 0; comp < map2d.numComp(); ++comp)
      {
        double val = iter(comp);          //<----
        assert(val == iter.value(comp));  //another way to get the value
        sum += val;
      }
      assert(iter(0) == *iter); //2 ways to get the first component value
      assert(iter(0) == iter.value());
      
    }
  }
  cout << end_timer() << "\t";
  assert(x_sum == sum);

  cout << endl;
}


int main(int argc, char** argv)
{
  //default options
  bool use_sparse = true;
  int nmats = 50;
  int ncells = 2000;
  int ncomp = 2;


  various_traversal_methods(nmats, ncells, ncomp, use_sparse);

  cout << "\nPress Enter to terminate...";
  cin.ignore();
}
