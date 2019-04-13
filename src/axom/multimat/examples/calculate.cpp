// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file calculate.cpp
 *
 * \brief Examples using MultiMat to do some calculation common in physics
 * simulation.
 */

#include "axom/multimat/multimat.hpp"


#include "axom/slic/core/UnitTestLogger.hpp"
#include "axom/core/utilities/Timer.hpp"

#include <vector>

#include "helper.hpp"
#include <ctime>

using namespace std;
using namespace axom::multimat;

#define ITERMAX 10  //define how many iterations to run test code

float filled_fraction = -1;

axom::utilities::Timer timer;

struct Robey_data
{
  /*
  * Data structure to store Robey's data.
  * When constructed, the volume fraction data is created. 
  * Then a layout is selected for the other array to be filled in.
  */

  int ncells;
  int nmats;
  float filled_percentage;

  //For creating MultiMat object. Always in dense cell-dominant layout.
  vector<bool> Volfrac_bool;
  vector<double> Volfrac_CD; //cell-dominant full volfrac array


  vector<double> Vol; //per cell, for all layouts

  //Per cellmat
  vector<double> Volfrac;
  vector<double> Densityfrac;
  vector<double> Temperaturefrac;
  vector<double> Pressurefrac;


  int cellmatcount;

  vector<double> nmatconsts; 

  int nnbrs_max;      //max number of neighbor = 8 for a 2d structured mesh
  vector<int> nnbrs;  //number of neighbors
  vector<int> nbrs;   //neighbor element id
  vector<double> cen; //centroids of cells

  //For CSR layout
  vector<int> begin_idx;
  vector<int> col_idx;
  vector<double> Volfrac_sparse; 
  vector<double> Densityfrac_sparse;
  vector<double> Temperaturefrac_sparse;
  vector<double> Pressurefrac_sparse;


  Robey_data(std::string filename = "", int ncells_in = 100, int nmats_in = 50)
  {

    if (filename != "")
    {
      //read from file... large and takes a long time.
      read_vol_frac_matrix_file(filename, ncells, nmats, Volfrac_CD,
                                filled_percentage);
    }
    else
    {
      //create random data
      //get_vol_frac_matrix_rand(ncells, nmats, Volfrac, filled_percentage);
      get_vol_frac_matrix_rand(ncells, nmats, Volfrac_CD, filled_percentage,
        ncells_in, nmats_in); //small version
    }

    filled_fraction = filled_percentage / 100.0f;

    // Some variables on neighbors
    //float L_f = read_from_file_bool ? 0.5 : 1.0;
    //// ave frac of nbrs containing material
    //int nnbrs_ave = 8;  // nearly so; 4000 boundary cells in 1 million cells
    //                    // in 3D, nnbrs_ave would be 26
    nnbrs_max = 8;
    // Build up list of neighbors for each cell
    // Assuming a 2D structured mesh, each cell will have a maximum of 8 nbrs
    nnbrs.resize(ncells);
    nbrs.resize(ncells*nnbrs_max);
    cen.resize(ncells * 2);

    nmatconsts.resize(nmats, 5.0);

    get_neighbors(ncells, nnbrs_max, nnbrs, nbrs);

    // Compute centroids of cells
    get_centroids(ncells, cen);

    //Making data for SLAM cell to mat relation
    Volfrac_bool.resize(ncells*nmats, false);
    cellmatcount = 0;
    for (unsigned int i = 0 ; i < Volfrac_CD.size() ; i++)
    {
      if (Volfrac_CD[i] > 0)
      {
        Volfrac_bool[i] = true;
        cellmatcount++;
      }
    }

  } //end constructor

  void set_up_cell_dom_data() {
    make_other_field_data_celldom( 
      ncells, nmats, Volfrac_CD, Volfrac, Vol, Densityfrac,
      Temperaturefrac, Pressurefrac, 
      Volfrac_sparse, Densityfrac_sparse, Temperaturefrac_sparse, Pressurefrac_sparse,
      begin_idx, col_idx);
  }
  void set_up_mat_dom_data() {
    make_other_field_data_matdom(
      ncells, nmats, Volfrac_CD, Volfrac, Vol, Densityfrac,
      Temperaturefrac, Pressurefrac, 
      Volfrac_sparse, Densityfrac_sparse, Temperaturefrac_sparse, Pressurefrac_sparse, 
      begin_idx, col_idx);
  }
};


//    Average density - Cell-Dominant Full Matrix
//    Robey's
void average_density_cell_dom_full(Robey_data& data){
  int ncells = data.ncells;
  int nmats = data.nmats;
  vector<double>& Volfrac = data.Volfrac;
  vector<double>& Densityfrac = data.Densityfrac;
  vector<double>& Vol = data.Vol;

  SLIC_INFO("-- Averaging Density cell-dominant array-access --");
  vector<double> Density_average(ncells);

  timer.reset();

  for (int iter = 0 ; iter < ITERMAX ; iter++)
  {
    timer.start();

    for (int ic = 0 ; ic < ncells ; ic++)
    {
      double density_ave = 0.0;
      for (int m = 0 ; m < nmats ; m++)
      {
        density_ave += Densityfrac[ic*nmats + m] * Volfrac[ic*nmats + m];
      }
      Density_average[ic] = density_ave / Vol[ic];
    }

    timer.stop();
  }

  double act_perf = timer.elapsed() / ITERMAX;
  SLIC_INFO("Average Density of mixed material cells    compute time is "
            << act_perf << " secs\n");

}

//    Average density - Cell-Dominant Compact
//    Robey's
void average_density_cell_dom_compact(Robey_data& data) {
  int ncells = data.ncells;
  int nmats = data.nmats;
  vector<double>& Volfrac = data.Volfrac_sparse;
  vector<double>& Densityfrac = data.Densityfrac_sparse;
  vector<double>& Vol = data.Vol;
  vector<int>& begin_idx = data.begin_idx;

  SLIC_INFO("-- Averaging Density cell-dominant compact array-access --");
  vector<double> Density_average(ncells);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; iter++)
  {
    timer.start();

    for (int ic = 0; ic < ncells; ic++)
    {
      double density_ave = 0.0;
      for (int ii = begin_idx[ic]; ii < begin_idx[ic + 1]; ii++)
      {
        density_ave += Densityfrac[ii] * Volfrac[ii];
      }
      Density_average[ic] = density_ave / Vol[ic];
    }

    timer.stop();
  }

  double act_perf = timer.elapsed() / ITERMAX;
  SLIC_INFO("Average Density of mixed material cells    compute time is "
    << act_perf << " secs\n");
}


//    Average density - Cell-Dominant Full Matrix
//    MultiMat - Direct Access
void average_density_cell_dom_mm_direct(MultiMat& mm) {
  SLIC_INFO("-- Averaging Density cell-dominant using MultiMat Direct (Dense) Access --");

  mm.convertLayoutToCellDominant();
  SLIC_INFO("MultiMat layout: " << mm.getDataLayoutAsString() << " & "
    << mm.getSparsityLayoutAsString());
  SLIC_ASSERT(mm.isDense());
  SLIC_ASSERT(mm.isCellDom());

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double> & Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double> & Volfrac = mm.get2dField<double>("Volfrac");
  MultiMat::Field1D<double> & Vol = mm.get1dField<double>("Vol");

  vector<double> Density_average(ncells);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; iter++)
  {
    timer.start();

    for (int ic = 0; ic < ncells; ic++)
    {
      double density_ave = 0.0;
      for (int m = 0; m < nmats; m++)
      {
        density_ave += *Densityfrac.findValue(ic, m) * *Volfrac.findValue(ic, m);
      }
      Density_average[ic] = density_ave / Vol(ic);
    }

    timer.stop();
  }
  double act_perf = timer.elapsed() / ITERMAX;
  SLIC_INFO("Average Density of mixed material cells    compute time is "
    << act_perf << " secs\n");
}

//    Average density - Cell-Dominant
//    MultiMat - Submap
void average_density_cell_dom_mm_submap(MultiMat& mm) {
  SLIC_INFO("-- Averaging Density cell-dominant using MultiMat Submap --");

  mm.convertLayoutToCellDominant();
  SLIC_INFO("MultiMat layout: " << mm.getDataLayoutAsString() << " & "
    << mm.getSparsityLayoutAsString());
  SLIC_ASSERT(mm.isCellDom());

  int ncells = mm.getNumberOfCells();
  MultiMat::Field2D<double> & Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double> & Volfrac = mm.get2dField<double>("Volfrac");
  MultiMat::Field1D<double> & Vol = mm.get1dField<double>("Vol");

  vector<double> Density_average(ncells);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; iter++)
  {
    timer.start();

    for (int ic = 0; ic < ncells; ic++)
    {
      double density_ave = 0.0;
      MultiMat::SubField<double> Densityfrac_row = Densityfrac(ic);
      MultiMat::SubField<double> Volfrac_row = Volfrac(ic);
      for (int j = 0; j < Densityfrac_row.size(); ++j)
      {
        density_ave += Densityfrac_row(j) * Volfrac_row(j);
      }
      Density_average[ic] = density_ave / Vol(ic);
    }

    timer.stop();
  }
  double act_perf = timer.elapsed() / ITERMAX;
  SLIC_INFO("Average Density of mixed material cells    compute time is "
    << act_perf << " secs\n");
}

//    Average density - Cell-Dominant
//    MultiMat - Index Array
void average_density_cell_dom_mm_idxarray(MultiMat& mm) {
  SLIC_INFO("-- Averaging Density cell-dominant using MultiMat Index Array --");

  mm.convertLayoutToCellDominant();
  SLIC_INFO("MultiMat layout: " << mm.getDataLayoutAsString() << " & "
    << mm.getSparsityLayoutAsString());
  SLIC_ASSERT(mm.isCellDom());

  int ncells = mm.getNumberOfCells();
  MultiMat::Field2D<double> & Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double> & Volfrac = mm.get2dField<double>("Volfrac");
  MultiMat::Field1D<double> & Vol = mm.get1dField<double>("Vol");

  vector<double> Density_average(ncells);
  
  timer.reset();

  for (int iter = 0; iter < ITERMAX; iter++)
  {
    timer.start();

    for (int ic = 0; ic < ncells; ic++)
    {
      auto idxSet = mm.getIndexingSetOfCell(ic);
      auto matId = mm.getMatInCell(ic);
      double density_ave = 0.0;
      for (int j = 0; j < idxSet.size(); ++j)
      {
        density_ave += Densityfrac[idxSet[j]] * Volfrac[idxSet[j]];
      }
      Density_average[ic] = density_ave / Vol(ic);
    }

    timer.stop();
  }
  double act_perf = timer.elapsed() / ITERMAX;
  SLIC_INFO("Average Density of mixed material cells    compute time is "
    << act_perf << " secs\n");
}

//    Average density - Cell-Dominant 
//    MultiMat - Flat Iterator
void average_density_mm_flatiter(MultiMat& mm) {
  SLIC_INFO("-- Averaging Density cell-dominant using MultiMat Flat Iter --");

  SLIC_INFO("MultiMat layout: " << mm.getDataLayoutAsString() << " & "
                                << mm.getSparsityLayoutAsString());
  //SLIC_ASSERT(mm.isCellDom());

  int ncells = mm.getNumberOfCells();
  MultiMat::Field2D<double>& Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double>& Volfrac = mm.get2dField<double>("Volfrac");
  MultiMat::Field1D<double>& Vol = mm.get1dField<double>("Vol");

  vector<double> Density_average(ncells, 0.0);

  timer.reset();

  for (int iter = 0 ; iter < ITERMAX ; iter++)
  {
    for (auto& v : Density_average)
      v = 0.0;

    timer.start();

    auto DensityIter = Densityfrac.begin();
    auto VolfracIter = Volfrac.begin();

    auto DensityIterEnd = Densityfrac.end();
    for ( ; DensityIter != DensityIterEnd ; DensityIter++, VolfracIter++)
    {
      Density_average[DensityIter.firstIndex()] += *DensityIter * *VolfracIter;
    }

    for (int ic = 0 ; ic < ncells ; ic++)
    {
      Density_average[ic] /= Vol(ic);
    }

    timer.stop();
  }
  double act_perf = timer.elapsed() / ITERMAX;
  SLIC_INFO("Average Density of mixed material cells    compute time is "
            << act_perf <<" secs\n");
}

//    Average density - Cell-Dominant 
//    MultiMat - Per-row Iterator
void average_density_cell_dom_mm_iter(MultiMat& mm) {
  SLIC_INFO("-- Averaging Density cell-dominant using MultiMat Iterator Per-row --");

  mm.convertLayoutToCellDominant();
  SLIC_INFO("MultiMat layout: " << mm.getDataLayoutAsString() << " & "
    << mm.getSparsityLayoutAsString());
  SLIC_ASSERT(mm.isCellDom());

  int ncells = mm.getNumberOfCells();
  MultiMat::Field2D<double> & Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double> & Volfrac = mm.get2dField<double>("Volfrac");
  MultiMat::Field1D<double> & Vol = mm.get1dField<double>("Vol");
  
  vector<double> Density_average(ncells, 0.0);
  
  timer.reset();

  for (int iter = 0; iter < ITERMAX; iter++)
  {
    for (auto& v : Density_average)
      v = 0.0;

    timer.start();
    
    for (int ic = 0; ic < ncells; ic++)
    {
      double density_ave = 0.0;
      auto DensityIter = Densityfrac.begin(ic);
      auto VolfracIter = Volfrac.begin(ic);
      auto DensityIterEnd = Densityfrac.end(ic);
      while (DensityIter != DensityIterEnd)
      {
        density_ave += *DensityIter * *VolfracIter;
        ++DensityIter; ++VolfracIter;
      }
      Density_average[ic] = density_ave / Vol(ic);
    }

    timer.stop();
  }
  double act_perf = timer.elapsed() / ITERMAX;
  SLIC_INFO("Average Density of mixed material cells    compute time is "
    << act_perf << " secs\n");
}



////////////////////// Material Dominant /////////////////////////

//    Average density - Material-Dominant Full Matrix
//    Robey's
double average_density_mat_dom_full(Robey_data& data) {
	int ncells = data.ncells;
	int nmats = data.nmats;
	vector<double>& Volfrac = data.Volfrac;
	vector<double>& Densityfrac = data.Densityfrac;
	vector<double>& Vol = data.Vol;

	SLIC_INFO("-- Averaging Density material-dominant full matrix array-access --");

  vector<double> Density_average(ncells, 0.0);
	
  timer.reset();

	for (int iter = 0; iter< ITERMAX; iter++) 
  {
    for (auto& v : Density_average)
      v = 0.0;

    timer.start();

		for (int m = 0; m < nmats; m++) {
			for (int ic = 0; ic < ncells; ic++) {
				Density_average[ic] += Densityfrac[m*ncells + ic] * Volfrac[m*ncells + ic];
			}
		}
		for (int ic = 0; ic < ncells; ic++) {
			Density_average[ic] /= Vol[ic];
		}

		timer.stop();
	}

	double act_perf = timer.elapsed() / ITERMAX;
	SLIC_INFO("Average Density of mixed material cells    compute time is "
		<< act_perf << " secs\n");
	return act_perf;
}

//    Average density - Material-Dominant Compact
//    Robey's
double average_density_mat_dom_compact(Robey_data& data) {
  int ncells = data.ncells;
  int nmats = data.nmats;
  vector<double>& Volfrac = data.Volfrac;
  vector<double>& Densityfrac = data.Densityfrac;
  vector<double>& Vol = data.Vol;
  vector<int>& begin_idx = data.begin_idx;
  vector<int>& cell_id = data.col_idx;

  SLIC_INFO("-- Averaging Density material-dominant compact array-access --");

  vector<double> Density_average(ncells, 0.0);
  
  timer.reset();

  for (int iter = 0; iter < ITERMAX; iter++) 
  {
    for (auto &v : Density_average)
      v = 0.0;
    
    timer.start();

    for (int m = 0; m < nmats; m++) {
      for (int ii = begin_idx[m]; ii < begin_idx[m + 1]; ii++)
      {
        Density_average[cell_id[ii]] += Densityfrac[ii] * Volfrac[ii];
      }
    }
    for (int ic = 0; ic < ncells; ic++) {
      Density_average[ic] /= Vol[ic];
    }

    timer.stop();
  }

  double act_perf = timer.elapsed() / ITERMAX;
  SLIC_INFO("Average Density of mixed material cells    compute time is "
    << act_perf << " secs\n");
  return act_perf;
}


//    Average density - Material-Dominant
//    MultiMat - Direct Access
void average_density_mat_dom_mm_direct(MultiMat& mm) {
  SLIC_INFO("-- Averaging Density mat-dominant using MultiMat --");

  mm.convertLayoutToMaterialDominant();

  SLIC_INFO("MultiMat layout: " << mm.getDataLayoutAsString() << " & "
    << mm.getSparsityLayoutAsString());
  SLIC_ASSERT(mm.isMatDom());
  SLIC_ASSERT(mm.isDense());

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double> & Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double> & Volfrac = mm.get2dField<double>("Volfrac");
  MultiMat::Field1D<double> & Vol = mm.get1dField<double>("Vol");

  vector<double> Density_average(ncells, 0.0);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; iter++)
  {
    for (auto& v : Density_average)
      v = 0.0;

    timer.start();

    for (int m = 0; m < nmats; ++m)
    {
      for (int ic = 0; ic < ncells; ++ic)
      {
        Density_average[ic] += *Densityfrac.findValue(m,ic) * *Volfrac.findValue(m,ic);
      }
    }
    for (int ic = 0; ic < ncells; ic++)
    {
      Density_average[ic] /= Vol(ic);
    }

    timer.stop();
  }
  double act_perf = timer.elapsed() / ITERMAX;

  SLIC_INFO("Average Density of mixed material cells    compute time is "
    << act_perf << " secs\n");
}

//    Average density - Material-Dominant
//    MultiMat - Submap
void average_density_mat_dom_mm_submap(MultiMat& mm) {
  SLIC_INFO("-- Averaging Density mat-dominant using MultiMat --");

  mm.convertLayoutToMaterialDominant();

  SLIC_INFO("MultiMat layout: " << mm.getDataLayoutAsString() << " & "
                                << mm.getSparsityLayoutAsString());
  SLIC_ASSERT(mm.isMatDom());

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double>& Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double>& Volfrac = mm.get2dField<double>("Volfrac");
  MultiMat::Field1D<double>& Vol = mm.get1dField<double>("Vol");

  vector<double> Density_average(ncells, 0.0);

  timer.reset();

  for (int iter = 0 ; iter < ITERMAX ; iter++)
  {
    for (auto& v : Density_average)
      v = 0.0;

    timer.start();

    for (int m = 0 ; m < nmats ; m++)
    {
      MultiMat::SubField<double> Densityfrac_row = Densityfrac(m);
      MultiMat::SubField<double> Volfrac_row = Volfrac(m);
      for (int j = 0 ; j < Volfrac_row.size() ; j++)
      {
        Density_average[j] += Densityfrac_row(j) * Volfrac_row(j);
      }
    }
    for (int ic = 0 ; ic < ncells ; ic++)
    {
      Density_average[ic] /= Vol(ic);
    }

    timer.stop();
  }
  double act_perf = timer.elapsed() / ITERMAX;

  SLIC_INFO("Average Density of mixed material cells    compute time is "
            << act_perf << " secs\n");
}

//    Average density - Material-Dominant
//    MultiMat - IndexArray
void average_density_mat_dom_mm_idxarray(MultiMat& mm) {
  SLIC_INFO("-- Averaging Density mat-dominant using MultiMat --");

  mm.convertLayoutToMaterialDominant();

  SLIC_INFO("MultiMat layout: " << mm.getDataLayoutAsString() << " & "
    << mm.getSparsityLayoutAsString());
  SLIC_ASSERT(mm.isMatDom());

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double> & Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double> & Volfrac = mm.get2dField<double>("Volfrac");
  MultiMat::Field1D<double> & Vol = mm.get1dField<double>("Vol");

  vector<double> Density_average(ncells, 0.0);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; iter++)
  {
    for (auto& v : Density_average)
      v = 0.0;

    timer.start();

    for (int m = 0; m < nmats; m++)
    {
      auto idxSet = mm.getIndexingSetOfMat(m);
      auto cellId = mm.getCellContainingMat(m);
      for (int j = 0; j < idxSet.size(); j++)
      {
        Density_average[cellId[j]] += Densityfrac[idxSet[j]] * Volfrac[idxSet[j]];
      }
    }
    for (int ic = 0; ic < ncells; ic++)
    {
      Density_average[ic] /= Vol(ic);
    }

    timer.stop();
  }
  double act_perf = timer.elapsed() / ITERMAX;

  SLIC_INFO("Average Density of mixed material cells    compute time is "
    << act_perf << " secs\n");
}

//    Average density - Material-Dominant
//    MultiMat - Iterator
void average_density_mat_dom_mm_iter(MultiMat& mm) {
  SLIC_INFO("-- Averaging Density mat-dominant using MultiMat --");

  mm.convertLayoutToMaterialDominant();

  SLIC_INFO("MultiMat layout: " << mm.getDataLayoutAsString() << " & "
    << mm.getSparsityLayoutAsString());
  SLIC_ASSERT(mm.isMatDom());

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double> & Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double> & Volfrac = mm.get2dField<double>("Volfrac");
  MultiMat::Field1D<double> & Vol = mm.get1dField<double>("Vol");

  vector<double> Density_average(ncells, 0.0);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; iter++)
  {
    for (auto& v : Density_average)
      v = 0.0;

    timer.start();

    for (int m = 0; m < nmats; m++)
    {
      auto Density_iter = Densityfrac.begin(m);
      auto Volfrac_iter = Volfrac.begin(m);
      auto Density_iterend = Densityfrac.end(m);
      while (Density_iter != Density_iterend)
      {
        Density_average[Density_iter.index()] += *Density_iter * *Volfrac_iter;
        ++Density_iter; ++Volfrac_iter;
      }
    }
    for (int ic = 0; ic < ncells; ic++)
    {
      Density_average[ic] /= Vol(ic);
    }

    timer.stop();
  }
  double act_perf = timer.elapsed() / ITERMAX;

  SLIC_INFO("Average Density of mixed material cells    compute time is "
    << act_perf << " secs\n");
}



//    Average density with if - Cell-Dominant Full Matrix
//      Robey's 
void average_density_cell_dom_with_if(Robey_data& data)
{
  int ncells = data.ncells;
  int nmats = data.nmats;
  vector<double>& Volfrac = data.Volfrac;
  vector<double>& Densityfrac = data.Densityfrac;
  vector<double>& Vol = data.Vol;

  SLIC_INFO("-- Averaging Density with if --");

  vector<double> Density_average(ncells, 0.0);
  timer.reset();

  for (int iter = 0 ; iter < ITERMAX ; iter++)
  {
    timer.start();

    for (int ic = 0 ; ic < ncells ; ic++)
    {
      double density_ave = 0.0;
      for (int m = 0 ; m < nmats ; m++)
      {
        if (Volfrac[ic*nmats + m] > 0.0)
        {
          density_ave += Densityfrac[ic*nmats + m] * Volfrac[ic*nmats + m];
        }
      }
      Density_average[ic] = density_ave / Vol[ic];
    }

    timer.stop();
  }

  double act_perf = timer.elapsed() / ITERMAX;
  SLIC_INFO("Average Density of frac with if            compute time is "
            << act_perf << " secs\n");

}

//    Average density with if - Cell-Dominant
//      MultiMat 
void average_density_cell_dom_with_if_mm(MultiMat&) {
  //SLIC_INFO("-- Averaging Density with if using MultiMat --");
  //this doesn't really make sense for a MultiMat stored in sparse layout...
}


//   Calculate pressure using ideal gas law - Cell-Dominant Full Matrix 
//     Robey's
void calculate_pressure(Robey_data& data)
{
  int ncells = data.ncells;
  int nmats = data.nmats;
  vector<double>& Volfrac = data.Volfrac;
  vector<double>& Densityfrac = data.Densityfrac;
  vector<double>& Temperaturefrac = data.Temperaturefrac;
  vector<double>& nmatconsts = data.nmatconsts;

  SLIC_INFO("-- Calculating pressure --");
  vector<double> Pressurefrac(ncells*nmats, 0);

  timer.reset();

  for (int iter = 0 ; iter< ITERMAX ; iter++)
  {
    timer.start();

    for (int ic = 0 ; ic < ncells ; ic++)
    {
      for (int m = 0 ; m < nmats ; m++)
      {
        if (Volfrac[ic*nmats + m] > 0.)
        {
          Pressurefrac[ic*nmats + m] = (nmatconsts[m] *
                                        Densityfrac[ic*nmats + m] *
                                        Temperaturefrac[ic*nmats + m])
                                       / (Volfrac[ic*nmats + m]);
        }
        else
        {
          Pressurefrac[ic*nmats + m] = 0.0;
        }
      }
    }

    timer.stop();
  }

  double act_perf = timer.elapsed() / ITERMAX;
  SLIC_INFO("Pressure Calculation of mixed material cells with if "
            << "compute time is " << act_perf << " secs\n");
}


//   Calculate pressure using ideal gas law - Cell-Dominant
//     MultiMat - Submap
void calculate_pressure_mm(MultiMat& mm)
{
  SLIC_INFO("-- Calculating pressure, using MultiMat --");
  mm.convertLayoutToCellDominant();
  SLIC_INFO("MultiMat layout: " << mm.getDataLayoutAsString() << " & "
                                << mm.getSparsityLayoutAsString());
  SLIC_ASSERT(mm.isCellDom());

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double>& Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double>& Volfrac = mm.get2dField<double>("Volfrac");
  MultiMat::Field2D<double>& Temperaturefrac =
    mm.get2dField<double>("Tempfrac");
  MultiMat::Field1D<double>& nmatconsts = mm.get1dField<double>("nmatconsts");
  
  vector<double> Pressurefrac(ncells*nmats, 0.0);

  timer.reset();

  for (int iter = 0 ; iter< ITERMAX ; iter++)
  {
    timer.start();

    for (int ic = 0 ; ic < ncells ; ic++)
    {
      MultiMat::SubField<double> Densityfrac_row = Densityfrac(ic);
      MultiMat::SubField<double> Volfrac_row = Volfrac(ic);
      MultiMat::SubField<double> Tempfrac_row = Temperaturefrac(ic);

      for (int j = 0 ; j < Volfrac_row.size() ; j++)
      {
        int m = Volfrac_row.index(j);
        Pressurefrac[ic*nmats + m] = ( nmatconsts(m) * Densityfrac_row(j) *
                                       Tempfrac_row(j)) / Volfrac_row(j);
      }
    }

    timer.stop();
  }

  double act_perf = timer.elapsed() / ITERMAX;
  SLIC_INFO("Pressure Calculation of mixed material cells with if "
            << "compute time is "<< act_perf << " secs\n");
}


//    Average material density over neighborhood of each cell
//      Robey's - Cell-Dominant Full Matrix
void average_material_density_over_cell_nbr(Robey_data& data)
{
  int ncells = data.ncells;
  int nmats = data.nmats;
  std::vector<double>& Volfrac = data.Volfrac;
  std::vector<double>& Densityfrac = data.Densityfrac;
  const std::vector<double>& cen = data.cen;
  const std::vector<int>& nnbrs = data.nnbrs;
  const std::vector<int>& nbrs = data.nbrs;


  SLIC_INFO("-- Calculating avg material density over cell neighbor --");
  timer.reset();

  for (int iter = 0 ; iter< ITERMAX ; iter++)
  {
    vector<double> MatDensity_average(ncells*nmats, 0);

    timer.start();

    for (int ic = 0 ; ic < ncells ; ic++)
    {
      double xc[2];
      xc[0] = cen[ic * 2]; xc[1] = cen[ic * 2 + 1];
      int nn = nnbrs[ic];
      int cnbrs[8];
      double dsqr[8];

      for (int n = 0 ; n < nn ; n++)
        cnbrs[n] = nbrs[ic * 8 + n];

      for (int n = 0 ; n < nn ; n++)
      {
        dsqr[n] = 0.0;
        for (int d = 0 ; d < 2 ; d++) //1???
        {
          double ddist = (xc[d] - cen[cnbrs[n] * 2 + d]);
          dsqr[n] += ddist * ddist;
        }
      }

      for (int m = 0 ; m < nmats ; m++)
      {
        if (Volfrac[ic*nmats + m] > 0.0)
        {
          int nnm = 0;         // number of nbrs with this material
          for (int n = 0 ; n < nn ; n++)
          {
            int jc = cnbrs[n];
            if (Volfrac[jc*nmats + m] > 0.0)
            {
              MatDensity_average[ic*nmats + m] +=
                Densityfrac[ic*nmats + m] / dsqr[n];
              nnm++;
            }
          }
          if(nnm>0) //This is not in Robey's code, but added to prevent NAN
            MatDensity_average[ic*nmats + m] /= nnm;
        }
        else
        {
          MatDensity_average[ic*nmats + m] = 0.0;
        }
      }
    }

    timer.stop();
  }

  double act_perf = timer.elapsed() / ITERMAX;

  SLIC_INFO("Average Material Density            compute time is "
            << act_perf << " secs\n");

}


//    Average material density over neighborhood of each cell
//      MultiMat - 
void average_material_density_over_cell_nbr_mm(MultiMat& mm, Robey_data& data)
{
  SLIC_INFO("-- Calculating avg material density over cell neighbor,"
            << " Multimat version --");
  mm.convertLayoutToCellDominant();
  SLIC_INFO("MultiMat layout: " << mm.getDataLayoutAsString() << " & "
                                << mm.getSparsityLayoutAsString());
  SLIC_ASSERT(mm.isCellDom());

  const std::vector<double>& cen = data.cen;
  const std::vector<int>& nnbrs = data.nnbrs;
  const std::vector<int>& nbrs = data.nbrs;

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double>& Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double>& Volfrac = mm.get2dField<double>("Volfrac");

  timer.reset();

  for (int iter = 0 ; iter< ITERMAX ; iter++)
  {
    vector<double> MatDensity_average(ncells*nmats, 0);
    timer.start();

    for (int ic = 0 ; ic < ncells ; ic++)
    {
      MultiMat::SubField<double> Densityfrac_row = Densityfrac(ic);
      MultiMat::SubField<double> Volfrac_row = Volfrac(ic);

      double xc[2];
      xc[0] = cen[ic * 2]; xc[1] = cen[ic * 2 + 1];
      int nn = nnbrs[ic];
      int cnbrs[8];
      double dsqr[8];

      for (int n = 0 ; n < nn ; n++)
        cnbrs[n] = nbrs[ic * 8 + n];

      for (int n = 0 ; n < nn ; n++)
      {
        dsqr[n] = 0.0;
        for (int d = 0 ; d < 2 ; d++) //1???
        {
          double ddist = (xc[d] - cen[cnbrs[n] * 2 + d]);
          dsqr[n] += ddist * ddist;
        }
      }

      for (int k = 0 ; k < Volfrac_row.size() ; k++)
      {
        int m = Volfrac_row.index(k);
        if (Volfrac_row(k)> 0.0)  //this check is not needed in sparse layout.
        {
          int nnm = 0;         // number of nbrs with this material
          for (int n = 0 ; n < nn ; n++)
          {
            int jc = cnbrs[n];
            double* jcvf = Volfrac.findValue(jc, m);
            if (jcvf != nullptr && *jcvf > 0.0)
            {
              MatDensity_average[ic*nmats + m] +=  Densityfrac_row(k) / dsqr[n];
              nnm++;
            }
          }
          if (nnm>0) //This is not in Robey's code, but added to prevent NAN
            MatDensity_average[ic*nmats + m] /= nnm;
        }
        else
        {
          MatDensity_average[ic*nmats + m] = 0.0;
        }
      }
    }

    timer.stop();
  }

  double act_perf = timer.elapsed() / ITERMAX;

  SLIC_INFO("Average Material Density            compute time is "
            << act_perf << " secs\n");

}


int main(int argc, char** argv)
{
  axom::slic::UnitTestLogger logger;

  if (argc != 1 && argc != 2 && argc != 3)
  {
    SLIC_WARNING("Usage: ./multimat_calculate_ex [<volfrac_data>]");
    return 1;
  }

  std::string fileName = "";
  int ncells_to_gen = 400;
  int nmats_to_gen = 20;

  if (argc == 2)
  {
    fileName = std::string(argv[1]);
  }
  else if(argc == 3)
  {
    ncells_to_gen = atoi(argv[1]);
    nmats_to_gen = atoi(argv[2]);
  }

  Robey_data data(fileName, ncells_to_gen, nmats_to_gen);
  data.set_up_cell_dom_data();

  //Set-up the multimat class
  MultiMat mm(DataLayout::CELL_CENTRIC, SparsityLayout::DENSE);
  mm.setNumberOfMaterials(data.nmats);
  mm.setNumberOfCells(data.ncells);
  mm.setCellMatRel(data.Volfrac_bool);

  //Setting field data in terms of slam
  mm.addField<>("Densityfrac", FieldMapping::PER_CELL_MAT,
                &data.Densityfrac[0]);
  mm.addField<>("Vol", FieldMapping::PER_CELL, &data.Vol[0]);
  mm.addField<>("Volfrac", FieldMapping::PER_CELL_MAT, &data.Volfrac[0]);
  mm.addField<>("Tempfrac", FieldMapping::PER_CELL_MAT,
                &data.Temperaturefrac[0]);
  mm.addField<>("Pressurefrac", FieldMapping::PER_CELL_MAT,
                &data.Pressurefrac[0]);
  mm.addField<>("nmatconsts", FieldMapping::PER_MAT, &data.nmatconsts[0]);

  //printself and check
  mm.isValid(true);

  //Run the Full layout, cell dom
  SLIC_INFO("*************** Full Layout - Cell Dominant **************");
  average_density_cell_dom_full(data);
  average_density_cell_dom_mm_direct(mm);
  average_density_cell_dom_mm_submap(mm);
  average_density_cell_dom_mm_iter(mm);
  //average_density_cell_dom_mm_idxarray(mm);
  average_density_mm_flatiter(mm);
  
  /*
  average_density_cell_dom_with_if(data);
  average_density_cell_dom_with_if_mm(mm);

  calculate_pressure(data);
  calculate_pressure_mm(mm);

  average_material_density_over_cell_nbr(data);
  average_material_density_over_cell_nbr_mm(mm, data);
  */

  //Run the Compact layout, cell dom
  SLIC_INFO("*************** Compact Layout - Cell Dominant **************");
  average_density_cell_dom_compact(data);
  mm.convertLayoutToSparse();
  average_density_cell_dom_mm_submap(mm);
  average_density_cell_dom_mm_iter(mm);
  average_density_cell_dom_mm_idxarray(mm);
  average_density_mm_flatiter(mm);


  //Run the Full layout, material dom
  SLIC_INFO("*************** Dense Layout - Material Dominant **************");
  data.set_up_mat_dom_data();
  mm.convertLayoutToDense();
  mm.convertLayoutToMaterialDominant();

  average_density_mat_dom_full(data);

  average_density_mat_dom_mm_direct(mm);
  average_density_mat_dom_mm_submap(mm);
  //average_density_mat_dom_mm_idxarray(mm);
  average_density_mat_dom_mm_iter(mm);
  average_density_mm_flatiter(mm);
  
  SLIC_INFO("*************** Compact Layout - Material Dominant **************");
  mm.convertLayoutToSparse();

  average_density_mat_dom_compact(data); 
  average_density_mat_dom_mm_submap(mm);
  average_density_mat_dom_mm_idxarray(mm);
  average_density_mat_dom_mm_iter(mm);
  average_density_mm_flatiter(mm);

  SLIC_INFO("*  Finished.");
  SLIC_INFO("**********************************************************");
  SLIC_INFO("");

  return 0;
}
