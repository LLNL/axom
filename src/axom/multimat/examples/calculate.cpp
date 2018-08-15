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
 * \file calculate.cpp
 *
 * \brief Examples using MultiMat to do some calculation common in physics 
 * simulation. 
 */

#include "axom/multimat/multimat.hpp"


#include "axom/slic/core/UnitTestLogger.hpp"
#include "axom/core/utilities/Timer.hpp"

#include "helper.hpp"
#include <ctime>

using namespace std;
using namespace axom::multimat;

#define ITERMAX 10  //define how many iterations to run test code


int method = 0;     // VF initialization: 0 - random, 1 - read volfrac.dat
float filled_fraction = -1;

void various_traversal_methods();

axom::utilities::Timer timer;

struct Robey_data
{
  int ncells;
  int nmats;
  float filled_percentage;

  vector<bool> Volfrac_bool;

  vector<double> Volfrac;
  vector<double> Vol;
  vector<double> Densityfrac;
  vector<double> Temperaturefrac;
  vector<double> Pressurefrac;
  vector<double> nmatconsts;

  int cellmatcount;

  int nnbrs_max;      //max number of neighbor = 8 for a 2d structured mesh
  vector<int> nnbrs;  //number of neighbors
  vector<int> nbrs;   //neighbor element id
  vector<double> cen; //centroids of cells

  vector<double> Densityfrac_sparse;
  vector<double> Volfrac_sparse;
  vector<double> Temperaturefrac_sparse;
  vector<double> Pressurefrac_sparse;

  Robey_data(std::string filename = "")
  {

    if (filename != "") {
      //read from file... large and takes a long time.
      read_vol_frac_matrix_file(filename, ncells, nmats, Volfrac, filled_percentage);
    }
    else {
      //create random data
      //get_vol_frac_matrix_rand(ncells, nmats, Volfrac, filled_percentage); 
      get_vol_frac_matrix_rand(ncells, nmats, Volfrac, filled_percentage, 100, 50); //small version
    }
    make_other_field_data(ncells, nmats, Volfrac, Vol, Densityfrac, Temperaturefrac, Pressurefrac, nmatconsts);
    filled_fraction = filled_percentage / 100.0f;

    // Some variables on neighbors
    //float L_f = read_from_file_bool ? 0.5 : 1.0;  // ave frac of nbrs containing material
    //int nnbrs_ave = 8;  // nearly so; 4000 boundary cells in 1 million cells
                        //                  // in 3D, nnbrs_ave would be 26
    nnbrs_max = 8;
    // Build up list of neighbors for each cell
    // Assuming a 2D structured mesh, each cell will have a maximum of 8 nbrs
    nnbrs.resize(ncells);
    nbrs.resize(ncells*nnbrs_max);
    cen.resize(ncells * 2);

    get_neighbors(ncells, nnbrs_max, nnbrs, nbrs);

    // Compute centroids of cells
    get_centroids(ncells, cen);

    //Making data for SLAM cell to mat relation
    Volfrac_bool.resize(ncells*nmats, false);
    cellmatcount = 0;
    for (unsigned int i = 0; i < Volfrac.size(); i++) {
      if (Volfrac[i] > 0) {
        Volfrac_bool[i] = true;
        cellmatcount++;
      }
    }

    //Create the sparse version of the data
    Densityfrac_sparse.resize(cellmatcount);
    Volfrac_sparse.resize(cellmatcount);
    Temperaturefrac_sparse.resize(cellmatcount);
    Pressurefrac_sparse.resize(cellmatcount);

    int ii = 0;
    for (int ci = 0; ci < ncells; ci++) {
      for (int mi = 0; mi < nmats; mi++) {
        double v = Volfrac[ci*nmats + mi];
        if (v != 0) {
          Densityfrac_sparse[ii] = Densityfrac[ci*nmats + mi];
          Volfrac_sparse[ii] = Volfrac[ci*nmats + mi];
          Temperaturefrac_sparse[ii] = Temperaturefrac[ci*nmats + mi];
          Pressurefrac_sparse[ii] = Pressurefrac[ci*nmats + mi];
          ii++;
        }
      }
    }
  }
};



////////////////////////////////////////////////////////////////////////////////
//    Average density with fractional densities
//    Copied from Robey's code for Cell-Dominant Full Matrix Data Structure
////////////////////////////////////////////////////////////////////////////////
void average_density_cell_dom(Robey_data data){
  int ncells = data.ncells;
  int nmats = data.nmats;
  vector<double>& Volfrac = data.Volfrac;
  vector<double>& Densityfrac = data.Densityfrac;
  vector<double>& Vol = data.Vol;

  SLIC_INFO("-- Averaging Density cell-dominant array-access --");
  vector<double> Density_average(ncells);

  timer.reset();
  
  for (int iter = 0; iter < ITERMAX; iter++) {
    timer.start();

    for (int ic = 0; ic < ncells; ic++) {
      double density_ave = 0.0;
      for (int m = 0; m < nmats; m++) {
        density_ave += Densityfrac[ic*nmats + m] * Volfrac[ic*nmats + m];
      }
      Density_average[ic] = density_ave / Vol[ic];
    }

    timer.stop();
  }

  double act_perf = timer.elapsed();
  SLIC_INFO("Average Density of mixed material cells    compute time is "
    << act_perf << " secs\n");

}



////////////////////////////////////////////////////////////////////////////////
//    Average density with fractional densities
//    Same as the function above, but modified to use MultiMat class
////////////////////////////////////////////////////////////////////////////////
void average_density_cell_dom_mm(MultiMat& mm) {
  SLIC_INFO("-- Averaging Density cell-dominant using MultiMat --");
  
  mm.convertLayoutToCellDominant();
  SLIC_INFO("MultiMat layout: " << mm.getDataLayoutAsString() << " & " 
                                << mm.getSparcityLayoutAsString());
  SLIC_ASSERT(mm.isCellDom());

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double>& Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double>& Volfrac = mm.get2dField<double>("Volfrac");
  MultiMat::Field1D<double>& Vol = mm.get1dField<double>("Vol");

  vector<double> Density_average(ncells);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; iter++) {
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
  double act_perf = timer.elapsed();
  SLIC_INFO("Average Density of mixed material cells    compute time is "
    << act_perf <<" secs\n");
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    Average density with fractional densities
//        Material-first loop
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void average_density_mat_dom_mm(MultiMat& mm) {
  SLIC_INFO("-- Averaging Density mat-dominant using MultiMat --");

  mm.convertLayoutToMaterialDominant();

  SLIC_INFO("MultiMat layout: " << mm.getDataLayoutAsString() << " & "
    << mm.getSparcityLayoutAsString());
  SLIC_ASSERT(mm.isMatDom());

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double>& Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double>& Volfrac = mm.get2dField<double>("Volfrac");
  MultiMat::Field1D<double>& Vol = mm.get1dField<double>("Vol");

  timer.reset();

  for (int iter = 0; iter < ITERMAX; iter++) {
    vector<double> Density_average(ncells, 0.0);
    timer.start();
    
    for (int m = 0; m < nmats; m++) 
    {
      MultiMat::SubField<double> Densityfrac_row = Densityfrac(m);
      MultiMat::SubField<double> Volfrac_row = Volfrac(m);
      for (int j = 0; j < Volfrac_row.size(); j++)
      {
        Density_average[j] += Densityfrac_row(j) * Volfrac_row(j);
      }
    }
    for (int ic = 0; ic < ncells; ic++) {
      Density_average[ic] /= Vol(ic);
    }
    
    timer.stop();
  }
  double act_perf = timer.elapsed();

  SLIC_INFO("Average Density of mixed material cells    compute time is "
    << act_perf << " secs\n");
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    Average density with fractional densities with if test
//    Copied from Robey's code for Cell-Dominant Full Matrix Data Structure
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void average_density_cell_dom_with_if(Robey_data data)
{
  int ncells = data.ncells;
  int nmats = data.nmats;
  vector<double>& Volfrac = data.Volfrac;
  vector<double>& Densityfrac = data.Densityfrac;
  vector<double>& Vol = data.Vol;

  SLIC_INFO("-- Averaging Density with if --");

  vector<double> Density_average(ncells, 0.0);
  timer.reset();

  for (int iter = 0; iter < ITERMAX; iter++) {
    timer.start();

    for (int ic = 0; ic < ncells; ic++) 
    {
      double density_ave = 0.0;
      for (int m = 0; m < nmats; m++) 
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

  double act_perf = timer.elapsed();
  SLIC_INFO("Average Density of frac with if            compute time is "
    << act_perf << " secs\n");

}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    Average density with fractional densities with if test
//          Same as the function above, but modified to use MultiMat class
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void average_density_cell_dom_with_if_mm(MultiMat& mm) {
  //SLIC_INFO("-- Averaging Density with if using MultiMat --");
  //this doesn't really make sense for a MultiMat stored in sparse layout...  
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   Calculate pressure using ideal gas law
//    Copied from Robey's code for Cell-Dominant Full Matrix Data Structure
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calculate_pressure(Robey_data data)
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

  for (int iter = 0; iter< ITERMAX; iter++) 
  {
    timer.start();

    for (int ic = 0; ic < ncells; ic++) 
    {
      for (int m = 0; m < nmats; m++) 
      {
        if (Volfrac[ic*nmats + m] > 0.) 
        {
          Pressurefrac[ic*nmats + m] = (nmatconsts[m] * 
            Densityfrac[ic*nmats + m] * Temperaturefrac[ic*nmats + m])
            / (Volfrac[ic*nmats + m]);
        }
        else {
          Pressurefrac[ic*nmats + m] = 0.0;
        }
      }
    }

    timer.stop();
  }

  double act_perf = timer.elapsed();
  SLIC_INFO("Pressure Calculation of mixed material cells with if "
    << "compute time is " << act_perf << " msecs\n");
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   Calculate pressure using ideal gas law
//          Same as the function above, but modified to use MultiMat class
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calculate_pressure_mm(MultiMat& mm)
{
  SLIC_INFO("-- Calculating pressure, using MultiMat --");
  mm.convertLayoutToCellDominant();
  SLIC_INFO("MultiMat layout: " << mm.getDataLayoutAsString() << " & "
    << mm.getSparcityLayoutAsString());
  SLIC_ASSERT(mm.isCellDom());

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double>& Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double>& Volfrac = mm.get2dField<double>("Volfrac");
  MultiMat::Field2D<double>& Temperaturefrac = mm.get2dField<double>("Temperaturefrac");
  MultiMat::Field1D<double>& nmatconsts = mm.get1dField<double>("nmatconsts");
  
  
  vector<double> Pressurefrac(ncells*nmats, 0.0);

  timer.reset();

  for (int iter = 0; iter< ITERMAX; iter++) 
  {
    timer.start();

    for (int ic = 0; ic < ncells; ic++) 
    {
      MultiMat::SubField<double> Densityfrac_row = Densityfrac(ic);
      MultiMat::SubField<double> Volfrac_row = Volfrac(ic);
      MultiMat::SubField<double> Tempfrac_row = Temperaturefrac(ic);

      for (int j = 0; j < Volfrac_row.size(); j++)
      {
        int m = Volfrac_row.index(j);
        Pressurefrac[ic*nmats + m] = ( nmatconsts(m) * Densityfrac_row(j) *
          Tempfrac_row(j)) / Volfrac_row(j);
      }
    }

    timer.stop();
  }

  double act_perf = timer.elapsed();
  SLIC_INFO("Pressure Calculation of mixed material cells with if "
    << "compute time is "<< act_perf << " msecs\n");
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    Average material density over neighborhood of each cell
//    Copied from Robey's code for Cell-Dominant Full Matrix Data Structure
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void average_material_density_over_cell_nbr(Robey_data data)
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

  for (int iter = 0; iter< ITERMAX; iter++) {
    vector<double> MatDensity_average(ncells*nmats, 0);

    timer.start();

    for (int ic = 0; ic < ncells; ic++) {
      double xc[2];
      xc[0] = cen[ic * 2]; xc[1] = cen[ic * 2 + 1];
      int nn = nnbrs[ic];
      int cnbrs[8];
      double dsqr[8];

      for (int n = 0; n < nn; n++)
        cnbrs[n] = nbrs[ic * 8 + n];

      for (int n = 0; n < nn; n++) {
        dsqr[n] = 0.0;
        for (int d = 0; d < 2; d++) { //1???
          double ddist = (xc[d] - cen[cnbrs[n] * 2 + d]);
          dsqr[n] += ddist * ddist;
        }
      }

      for (int m = 0; m < nmats; m++) {
        if (Volfrac[ic*nmats + m] > 0.0) {
          int nnm = 0;         // number of nbrs with this material
          for (int n = 0; n < nn; n++) {
            int jc = cnbrs[n];
            if (Volfrac[jc*nmats + m] > 0.0) {
              MatDensity_average[ic*nmats + m] += Densityfrac[ic*nmats + m] / dsqr[n];
              nnm++;
            }
          }
          if(nnm>0) //This is not in Robey's code, but added to prevent NAN
            MatDensity_average[ic*nmats + m] /= nnm;
        }
        else {
          MatDensity_average[ic*nmats + m] = 0.0;
        }
      }
    }

    timer.stop();
  }

  double act_perf = timer.elapsed();

  SLIC_INFO("Average Material Density            compute time is "
    << act_perf << " msecs\n",);

}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    Average material density over neighborhood of each cell
//    Same as the function above, but modified to use MultiMat class
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void average_material_density_over_cell_nbr_mm(MultiMat& mm, Robey_data& data)
{
  SLIC_INFO("-- Calculating avg material density over cell neighbor,"
    << " Multimat version --");
  mm.convertLayoutToCellDominant();
  SLIC_INFO("MultiMat layout: " << mm.getDataLayoutAsString() << " & "
    << mm.getSparcityLayoutAsString());
  SLIC_ASSERT(mm.isCellDom());

  const std::vector<double>& cen = data.cen;
  const std::vector<int>& nnbrs = data.nnbrs;
  const std::vector<int>& nbrs = data.nbrs;

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double>& Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double>& Volfrac = mm.get2dField<double>("Volfrac");
  
  timer.reset();

  for (int iter = 0; iter< ITERMAX; iter++) {
    vector<double> MatDensity_average(ncells*nmats, 0);
    timer.start();

    for (int ic = 0; ic < ncells; ic++) {
      MultiMat::SubField<double> Densityfrac_row = Densityfrac(ic);
      MultiMat::SubField<double> Volfrac_row = Volfrac(ic);

      double xc[2];
      xc[0] = cen[ic * 2]; xc[1] = cen[ic * 2 + 1];
      int nn = nnbrs[ic];
      int cnbrs[8];
      double dsqr[8];

      for (int n = 0; n < nn; n++)
        cnbrs[n] = nbrs[ic * 8 + n];

      for (int n = 0; n < nn; n++) {
        dsqr[n] = 0.0;
        for (int d = 0; d < 2; d++) { //1???
          double ddist = (xc[d] - cen[cnbrs[n] * 2 + d]);
          dsqr[n] += ddist * ddist;
        }
      }

      for (int k = 0; k < Volfrac_row.size(); k++) 
      {
        int m = Volfrac_row.index(k);
        if (Volfrac_row(k)> 0.0)  //this check is not needed in sparse layout.
        {
          int nnm = 0;         // number of nbrs with this material
          for (int n = 0; n < nn; n++) {
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
        else {
          MatDensity_average[ic*nmats + m] = 0.0;
        }
      }
    }

    timer.stop();
  }

  double act_perf = timer.elapsed();  // divide by (double)ITERMAX ? 

  SLIC_INFO("Average Material Density            compute time is "
    << act_perf << " msecs\n");

}


int main(int argc, char** argv)
{
  axom::slic::UnitTestLogger logger;

  if (argc != 1 && argc != 2)
  {
    SLIC_WARNING("Usage: ./multimat_calculate_ex [<volfrac_data>]");
    return 1;
  }
  
  std::string fileName = "";
  if (argc == 2)
  {
    fileName = std::string(argv[1]);
  }

  Robey_data data(fileName);

  //Set-up the multimat class
  MultiMat mm;
  mm.setNumberOfMat(data.nmats);
  mm.setNumberOfCell(data.ncells);
  mm.setCellMatRel(data.Volfrac_bool);

  //Setting field data in terms of slam
  mm.addField<>("Densityfrac"    , FieldMapping::PER_CELL_MAT, &data.Densityfrac_sparse[0]);
  mm.addField<>("Vol"            , FieldMapping::PER_CELL    , &data.Vol[0]);
  mm.addField<>("Volfrac"        , FieldMapping::PER_CELL_MAT, &data.Volfrac_sparse[0]);
  mm.addField<>("Temperaturefrac", FieldMapping::PER_CELL_MAT, &data.Temperaturefrac_sparse[0]);
  mm.addField<>("Pressurefrac"   , FieldMapping::PER_CELL_MAT, &data.Pressurefrac_sparse[0]);
  mm.addField<>("nmatconsts"     , FieldMapping::PER_MAT     , &data.nmatconsts[0]);

  //printself and check
  
  mm.isValid(true);

  //Run the examples
  average_density_cell_dom(data);
  average_density_cell_dom_mm(mm);
  average_density_mat_dom_mm(mm);

  average_density_cell_dom_with_if(data);
  average_density_cell_dom_with_if_mm(mm);

  calculate_pressure(data);
  calculate_pressure_mm(mm);

  average_material_density_over_cell_nbr(data);
  average_material_density_over_cell_nbr_mm(mm, data);

  SLIC_INFO("*  Finished.");

  return 0;
}
