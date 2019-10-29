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

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"

#ifdef AXOM_DEBUG
# define ITERMAX 1   //define how many iterations to run test code
#else
# define ITERMAX 20  //define how many iterations to run test code
#endif

#include "helper.hpp"

namespace slam = axom::slam;
using namespace axom::multimat;

#define run_slam_bivarmap

template<typename B>
using Field2DT = MultiMat::Field2D<double, B>;

template<DataLayout D, typename B>
using Field2DTempT = MultiMat::Field2DTemplated<double, D, B>;

template<typename B>
using BiVarMapT = slam::BivariateMap<double, B>;

template<typename FieldTypeT>
using GetFieldFunType = FieldTypeT(MultiMat:: *)(const std::string&);

template<typename FieldTypeT, typename BSet>
GetFieldFunType<FieldTypeT> getfieldfun()
{
  if constexpr (std::is_same< FieldTypeT, Field2DT<BSet>>::value) {
    SLIC_INFO("Using templated field2D");
    return &MultiMat::get2dField<double, BSet>;
  }
  else if constexpr (std::is_same< FieldTypeT, Field2DTempT<DataLayout::CELL_DOM, BSet>>::value)
  {
    SLIC_INFO("Using specialized templated-field2D");
    return &MultiMat::getTemplated2DField<double, DataLayout::CELL_DOM, BSet>;
  }
  else if constexpr (std::is_same< FieldTypeT, Field2DTempT<DataLayout::MAT_DOM, BSet>>::value)
  {
    SLIC_INFO("Using specialized templated-field2D");
    return &MultiMat::getTemplated2DField<double, DataLayout::MAT_DOM, BSet>;
  }
  else if constexpr (std::is_same< FieldTypeT, BiVarMapT<BSet> >::value)
  {
    SLIC_INFO("Using slam BivariateMap");
    return &MultiMat::get2dFieldAsSlamBivarMap<double, BSet>;
  }
}


multirun_timer timer;
Value_Checker data_checker;
Result_Store result_store;


////////////////////// Average density - Cell Dominant /////////////////////////

//    Average density - Cell-Dominant Full Matrix
//    Robey's
void average_density_cell_dom_full(Robey_data& data){
  int ncells = data.ncells;
  int nmats = data.nmats;
  std::vector<double>& Volfrac = data.Volfrac;
  std::vector<double>& Densityfrac = data.Densityfrac;
  std::vector<double>& Vol = data.Vol;

  SLIC_INFO("-- Averaging Density, Cell-Dominant Full Matrix Array Access --");
  std::vector<double> Density_average(ncells);

  timer.reset();

  for (int iter = 0 ; iter < ITERMAX ; ++iter)
  {
    timer.start();

    for (int ic = 0 ; ic < ncells ; ++ic)
    {
      double density_ave = 0.0;
      for (int m = 0 ; m < nmats ; ++m)
      {
        density_ave += Densityfrac[ic*nmats + m] * Volfrac[ic*nmats + m];
      }
      Density_average[ic] = density_ave / Vol[ic];
    }

    timer.record();
    data_checker.check(Density_average);
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::avg_density, DataLayout::CELL_DOM,
    SparsityLayout::DENSE, Result_Store::method_csr, act_perf);

  SLIC_INFO("Average Density                      compute time is "
            << act_perf << " secs\n");

}

//    Average density - Cell-Dominant Compact
//    Robey's
void average_density_cell_dom_compact(Robey_data& data) {
  int ncells = data.ncells;
  std::vector<double>& Volfrac = data.Volfrac_sparse;
  std::vector<double>& Densityfrac = data.Densityfrac_sparse;
  std::vector<double>& Vol = data.Vol;
  std::vector<int>& begin_idx = data.begin_idx;

  SLIC_INFO("-- Averaging Density cell-dominant compact array-access --");
  std::vector<double> Density_average(ncells);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    timer.start();

    for (int ic = 0; ic < ncells; ++ic)
    {
      double den = 0.0;
      for (int ii = begin_idx[ic]; ii < begin_idx[ic + 1]; ++ii)
      {
          den += Densityfrac[ii] * Volfrac[ii];
      }
      Density_average[ic] = den / Vol[ic];
    }

    timer.record();
    data_checker.check(Density_average);
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::avg_density, DataLayout::CELL_DOM,
    SparsityLayout::SPARSE, Result_Store::method_csr, act_perf);

  SLIC_INFO("Average Density                      compute time is "
    << act_perf << " secs\n");
}

//    Average density - Cell-Dominant Full Matrix
//    MultiMat - Direct Access
void average_density_cell_dom_mm_template(MultiMat& mm) {
  SLIC_INFO("-- Averaging Density cell-dominant using MultiMat Direct (Dense) Access w/ template cast --");

  mm.convertLayoutToCellDominant(); 
  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldSparsityLayout(0) == SparsityLayout::DENSE);
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::CELL_DOM);

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  auto& Densityfrac = mm.get2dField<double>("Densityfrac");
  auto& Volfrac = mm.get2dField<double>("Volfrac");
  auto& Vol = mm.get1dField<double>("Vol");

  using SetType = slam::RangeSet<>;
  using ProductSetType = slam::ProductSet<SetType,SetType>;
  auto prodset = Densityfrac.getBivariateSet<ProductSetType>();

  std::vector<double> Density_average(ncells);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    timer.start();

    for (int ic = 0; ic < ncells; ++ic)
    {
      double density_ave = 0.0;
      for (int m = 0; m < nmats; ++m)
      {
         auto flatIdx = prodset.findElementFlatIndex(ic, m);
         density_ave += Densityfrac[flatIdx] * Volfrac[flatIdx];
      }
      Density_average[ic] = density_ave / Vol[ic];
    }

    timer.record();
    data_checker.check(Density_average);
  }
  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::avg_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_direct, act_perf);
  SLIC_INFO("Average Density                      compute time is "
    << act_perf << " secs\n");
}

//    Average density - Cell-Dominant Full Matrix
//    MultiMat - Direct Access
void average_density_cell_dom_mm_direct(MultiMat& mm) {
  SLIC_INFO("-- Averaging Density cell-dominant using MultiMat Direct (Dense) Access --");

  mm.convertLayoutToCellDominant();
  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldSparsityLayout(0) == SparsityLayout::DENSE);
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::CELL_DOM);

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  auto& Densityfrac = mm.get2dField<double>("Densityfrac");
  auto& Volfrac = mm.get2dField<double>("Volfrac");
  auto& Vol = mm.get1dField<double>("Vol");

  std::vector<double> Density_average(ncells);

  using SetType = slam::RangeSet<>;
  using ProductSet = slam::ProductSet<SetType,SetType>;
  using Stride = slam::policies::StrideOne<int>;
                 //slam::policies::RuntimeStride<int> >
  using Ind = slam::policies::STLVectorIndirection<int, double>;

  /**** BEGIN MOVE THIS CODE INTO MULTIMAT  *****/
  // convert DensityFrac to use ProductSets
  using BMapT = slam::BivariateMap<double, ProductSet, Ind, Stride>;

  auto prodSet = Densityfrac.getBivariateSet<ProductSet>();
  auto densityMap = BMapT(&prodSet);
  densityMap.copy(Densityfrac.getMap()->data().data() );

  auto volfracMap = BMapT(&prodSet);
  volfracMap.copy(Volfrac.getMap()->data().data() );
  /**** END MOVE THIS CODE INTO MULTIMAT  *****/


  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    timer.start();

    for (int ic = 0; ic < ncells; ++ic)
    {
      double density_ave = 0.0;
      for (int m = 0; m < nmats; ++m)
      {
        density_ave += densityMap(ic, m) * volfracMap(ic, m);
      }
      Density_average[ic] = density_ave / Vol[ic];
    }

    timer.record();
    data_checker.check(Density_average);
  }
  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::avg_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_direct, act_perf);
  SLIC_INFO("Average Density                      compute time is "
    << act_perf << " secs\n");
}

//    Average density - Cell-Dominant
//    MultiMat - Submap
void average_density_cell_dom_mm_submap(MultiMat& mm) {
  SLIC_INFO("-- Averaging Density cell-dominant using MultiMat Submap --");

  mm.convertLayoutToCellDominant(); 
  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::CELL_DOM); 

  int ncells = mm.getNumberOfCells();
  auto& Densityfrac = mm.get2dField<double>("Densityfrac");
  auto& Volfrac = mm.get2dField<double>("Volfrac");
  auto& Vol = mm.get1dField<double>("Vol");

  std::vector<double> Density_average(ncells);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    timer.start();

    for (int ic = 0; ic < ncells; ++ic)
    {
      double density_ave = 0.0;
      auto Densityfrac_row = Densityfrac(ic);
      auto Volfrac_row = Volfrac(ic);

      const auto cellMats = Densityfrac_row.size();
      for (int j = 0; j < cellMats; ++j)
      {
        density_ave += Densityfrac_row(j) * Volfrac_row(j);
      }
      Density_average[ic] = density_ave / Vol[ic];
    }

    timer.record();
    data_checker.check(Density_average);
  }
  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::avg_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_submap, act_perf);
  SLIC_INFO("Average Density                      compute time is "
    << act_perf << " secs\n");
}

//    Average density - Cell-Dominant
//    MultiMat - Submap
template<typename BSet>
void average_density_cell_dom_mm_submap(MultiMat& mm)
{
  SLIC_INFO("-- Averaging Density cell-dominant using MultiMat Submap -- templated on concrete bivariate set");

  mm.convertLayoutToCellDominant();
  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::CELL_DOM);

  int ncells = mm.getNumberOfCells();


  SLIC_INFO("Using slam bivariate set");
  auto Densityfrac = mm.get2dFieldAsSlamBivarMap<double, BSet>("Densityfrac");
  auto Volfrac = mm.get2dFieldAsSlamBivarMap<double, BSet>("Volfrac");

  auto& Vol = mm.get1dField<double>("Vol");

  std::vector<double> Density_average(ncells);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    timer.start();

    for (int ic = 0; ic < ncells; ++ic)
    {
      double density_ave = 0.0;
      auto Densityfrac_row = Densityfrac(ic);
      auto Volfrac_row = Volfrac(ic);

      const auto cellMats = Densityfrac_row.size();
      for (int j = 0; j < cellMats; ++j)
      {
        density_ave += Densityfrac_row(j) * Volfrac_row(j);
      }
      Density_average[ic] = density_ave / Vol[ic];
    }

    timer.record();
    data_checker.check(Density_average);
  }
  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::avg_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_submap_templated, act_perf);
  SLIC_INFO("Average Density                      compute time is "
    << act_perf << " secs\n");
}


//    Average density - Cell-Dominant
//    MultiMat - Submap
template<typename BSet, typename FieldTypeT >
void average_density_cell_dom_mm_submap(MultiMat& mm)
{
  SLIC_INFO("-- Averaging Density cell-dominant using MultiMat Submap -- templated on concrete bivariate set");

  mm.convertLayoutToCellDominant();
  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
      << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::CELL_DOM);

  int ncells = mm.getNumberOfCells();

  //pick which get field function to use
  GetFieldFunType<FieldTypeT> getFieldFun = getfieldfun<FieldTypeT, BSet>();

  /*if constexpr (std::is_same< FieldTypeT, Field2DT<BSet>>::value) {
    SLIC_INFO("Using templated field2D");
    getFieldFun = &MultiMat::get2dField<double, BSet>;
  }
  else if constexpr (std::is_same< FieldTypeT, Field2DTempT<DataLayout::CELL_DOM, BSet>>::value)
  {
    SLIC_INFO("Using specialized templated-field2D");
    getFieldFun = &MultiMat::getTemplated2DField<double, DataLayout::CELL_DOM, BSet>;
  }
  else if constexpr (std::is_same< FieldTypeT, BiVarMapT<BSet> >::value)
  {
    SLIC_INFO("Using slam BivariateMap");
    getFieldFun = &MultiMat::get2dFieldAsSlamBivarMap<double, BSet>;
  }*/


  FieldTypeT Densityfrac = (mm.*getFieldFun)("Densityfrac");
  FieldTypeT Volfrac = (mm.*getFieldFun)("Volfrac");

  auto& Vol = mm.get1dField<double>("Vol");

  std::vector<double> Density_average(ncells);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    timer.start();

    for (int ic = 0; ic < ncells; ++ic)
    {
      double density_ave = 0.0;
      auto Densityfrac_row = Densityfrac(ic);
      auto Volfrac_row = Volfrac(ic);

      const auto cellMats = Densityfrac_row.size();
      for (int j = 0; j < cellMats; ++j)
      {
        density_ave += Densityfrac_row(j) * Volfrac_row(j);
      }
      Density_average[ic] = density_ave / Vol[ic];
   }

    timer.record();
    data_checker.check(Density_average);
  }
  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::avg_density, mm.getFieldDataLayout(0),
      mm.getFieldSparsityLayout(0), Result_Store::mm_submap_templated, act_perf);
  SLIC_INFO("Average Density                      compute time is "
      << act_perf << " secs\n");
}



//    Average density - Cell-Dominant
//    MultiMat - Index Array
void average_density_cell_dom_mm_idxarray(MultiMat& mm) {
  SLIC_INFO("-- Averaging Density cell-dominant using MultiMat Index Array --");

  mm.convertLayoutToCellDominant();
  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::CELL_DOM);

  int ncells = mm.getNumberOfCells();
  auto& Densityfrac = mm.get2dField<double>("Densityfrac");
  auto& Volfrac = mm.get2dField<double>("Volfrac");
  auto& Vol = mm.get1dField<double>("Vol");

  std::vector<double> Density_average(ncells);
  
  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    timer.start();

    for (int ic = 0; ic < ncells; ++ic)
    {
      auto idxSet = mm.getIndexingSetOfCell(ic, mm.getFieldSparsityLayout(0));
      //auto matId = mm.getMatInCell(ic);
      double density_ave = 0.0;
      int sz = idxSet.size();
      for (int j = 0; j < sz; ++j)
      {
        int jj = idxSet[j];
        density_ave += Densityfrac[jj] * Volfrac[jj];
      }
      Density_average[ic] = density_ave / Vol[ic];
    }

    timer.record();
    data_checker.check(Density_average);
  }
  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::avg_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_idxarray, act_perf);
  SLIC_INFO("Average Density                      compute time is "
    << act_perf << " secs\n");
}

//    Average density - Cell-Dominant 
//    MultiMat - Flat Iterator
void average_density_cell_dom_mm_flatiter(MultiMat& mm) {
  SLIC_INFO("-- Averaging Density cell-dominant using MultiMat Flat Iter --");

  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
                                << mm.getFieldSparsityLayoutAsString(0));
  
  int ncells = mm.getNumberOfCells();
  MultiMat::Field2D<double>& Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double>& Volfrac = mm.get2dField<double>("Volfrac");
  MultiMat::Field1D<double>& Vol = mm.get1dField<double>("Vol");

  std::vector<double> Density_average(ncells, 0.0);

  timer.reset();

  for (int iter = 0 ; iter < ITERMAX ; ++iter)
  {
    for (auto& v : Density_average)
      v = 0.0;

    timer.start();

    auto DensityIter = Densityfrac.begin();
    auto VolfracIter = Volfrac.begin();

    auto DensityIterEnd = Densityfrac.end();
    for ( ; DensityIter != DensityIterEnd ; ++DensityIter, ++VolfracIter)
    {
      Density_average[DensityIter.firstIndex()] += *DensityIter * *VolfracIter;
    }

    for (int ic = 0 ; ic < ncells ; ++ic)
    {
      Density_average[ic] /= Vol[ic];
    }

    timer.record();
    data_checker.check(Density_average);
  }
  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::avg_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_flatiter, act_perf);
  SLIC_INFO("Average Density                      compute time is "
            << act_perf <<" secs\n");
}

//    Average density - Cell-Dominant 
//    MultiMat - Per-row Iterator
void average_density_cell_dom_mm_iter(MultiMat& mm) {
  SLIC_INFO("-- Averaging Density cell-dominant using MultiMat Iterator Per-row --");

  mm.convertLayoutToCellDominant();
  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::CELL_DOM);

  int ncells = mm.getNumberOfCells();
  MultiMat::Field2D<double> & Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double> & Volfrac = mm.get2dField<double>("Volfrac");
  MultiMat::Field1D<double> & Vol = mm.get1dField<double>("Vol");
  
  std::vector<double> Density_average(ncells, 0.0);
  
  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    for (auto& v : Density_average)
      v = 0.0;

    timer.start();
    
    for (int ic = 0; ic < ncells; ++ic)
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
      Density_average[ic] = density_ave / Vol[ic];
    }

    timer.record();
    data_checker.check(Density_average);
  }
  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::avg_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_iter, act_perf);
  SLIC_INFO("Average Density                      compute time is "
    << act_perf << " secs\n");
}



////////////////////// Average density - Material Dominant /////////////////////////

//    Average density - Material-Dominant Full Matrix
//    Robey's
double average_density_mat_dom_full(Robey_data& data) 
{
  int ncells = data.ncells;
  int nmats = data.nmats;
  std::vector<double>& Volfrac = data.Volfrac;
  std::vector<double>& Densityfrac = data.Densityfrac;
  std::vector<double>& Vol = data.Vol;
  
  SLIC_INFO("-- Averaging Density material-dominant full matrix array-access --");

  std::vector<double> Density_average(ncells, 0.0);
	
  timer.reset();

	for (int iter = 0; iter< ITERMAX; ++iter) 
  {
    for (auto& v : Density_average)
      v = 0.0;

    timer.start();

		for (int m = 0; m < nmats; ++m) {
			for (int ic = 0; ic < ncells; ++ic) {
        Density_average[ic] += Densityfrac[m*ncells + ic] * Volfrac[m*ncells + ic];
			}
		}
		for (int ic = 0; ic < ncells; ++ic) {
			Density_average[ic] /= Vol[ic];
		}

		timer.record();
    data_checker.check(Density_average);
	}

	double act_perf = timer.get_median();
   result_store.add_result(Result_Store::avg_density, DataLayout::MAT_DOM,
     SparsityLayout::DENSE, Result_Store::method_csr, act_perf);
	SLIC_INFO("Average Density                      compute time is "
		<< act_perf << " secs\n");
	return act_perf;
}

//    Average density - Material-Dominant Compact
//    Robey's
double average_density_mat_dom_compact(Robey_data& data) {
  int ncells = data.ncells;
  int nmats = data.nmats;
  std::vector<double>& Volfrac = data.Volfrac_sparse;
  std::vector<double>& Densityfrac = data.Densityfrac_sparse;
  std::vector<double>& Vol = data.Vol;
  std::vector<int>& begin_idx = data.begin_idx;
  std::vector<int>& cell_id = data.col_idx;

  SLIC_INFO("-- Averaging Density material-dominant compact array-access --");

  std::vector<double> Density_average(ncells, 0.0);
  
  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter) 
  {
    for (auto &v : Density_average)
      v = 0.0;
    
    timer.start();

    for (int m = 0; m < nmats; ++m) {
      for (int ii = begin_idx[m]; ii < begin_idx[m + 1]; ++ii)
      {
        Density_average[cell_id[ii]] += Densityfrac[ii] * Volfrac[ii];
      }
    }
    for (int ic = 0; ic < ncells; ++ic) {
      Density_average[ic] /= Vol[ic];
    }

    timer.record();
    data_checker.check(Density_average);
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::avg_density, DataLayout::MAT_DOM,
    SparsityLayout::SPARSE, Result_Store::method_csr, act_perf);
  SLIC_INFO("Average Density                      compute time is "
    << act_perf << " secs\n");
  return act_perf;
}


//    Average density - Material-Dominant
//    MultiMat - Direct Access
void average_density_mat_dom_mm_direct(MultiMat& mm) {
  SLIC_INFO("-- Averaging Density mat-dominant using MultiMat Direct Access --");

  mm.convertLayoutToMaterialDominant();

  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::MAT_DOM);
  SLIC_ASSERT(mm.getFieldSparsityLayout(0) == SparsityLayout::DENSE); 

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double> & Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double> & Volfrac = mm.get2dField<double>("Volfrac");
  MultiMat::Field1D<double> & Vol = mm.get1dField<double>("Vol");

  std::vector<double> Density_average(ncells, 0.0);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
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
    for (int ic = 0; ic < ncells; ++ic)
    {
      Density_average[ic] /= Vol[ic];
    }

    timer.record();
    data_checker.check(Density_average);
  }
  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::avg_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_direct, act_perf);

  SLIC_INFO("Average Density                      compute time is "
    << act_perf << " secs\n");
}

//    Average density - Material-Dominant
//    MultiMat - Submap
void average_density_mat_dom_mm_submap(MultiMat& mm) {
  SLIC_INFO("-- Averaging Density mat-dominant using MultiMat Submap --");

  mm.convertLayoutToMaterialDominant();

  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
                                << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::MAT_DOM);

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  auto& Densityfrac = mm.get2dField<double>("Densityfrac");
  auto& Volfrac = mm.get2dField<double>("Volfrac");
  auto& Vol = mm.get1dField<double>("Vol");

  std::vector<double> Density_average(ncells, 0.0);

  timer.reset();

  for (int iter = 0 ; iter < ITERMAX ; ++iter)
  {
    for (auto& v : Density_average)
      v = 0.0;

    timer.start();

    for (int m = 0 ; m < nmats ; ++m)
    {
      auto Densityfrac_row = Densityfrac(m);
      auto Volfrac_row = Volfrac(m);
      for (int j = 0 ; j < Volfrac_row.size() ; ++j)
      {
        Density_average[Densityfrac_row.index(j)] += Densityfrac_row(j) * Volfrac_row(j);
      }
    }
    for (int ic = 0 ; ic < ncells ; ++ic)
    {
      Density_average[ic] /= Vol[ic];
    }

    timer.record();
    data_checker.check(Density_average);
  }
  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::avg_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_submap, act_perf);

  SLIC_INFO("Average Density                      compute time is "
            << act_perf << " secs\n");
}


//    Average density - Material-Dominant
//    MultiMat - Submap
template<typename BSet, typename FieldTypeT>
void average_density_mat_dom_mm_submap(MultiMat& mm) {
  SLIC_INFO("-- Averaging Density mat-dominant using MultiMat Submap -- templated on concrete bivariate set");

  mm.convertLayoutToMaterialDominant();

  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::MAT_DOM);

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();

  //pick which get field function to use
  GetFieldFunType<FieldTypeT> getFieldFun = getfieldfun<FieldTypeT, BSet>();

  FieldTypeT Densityfrac = (mm.*getFieldFun)("Densityfrac");
  FieldTypeT Volfrac = (mm.*getFieldFun)("Volfrac");

  //auto Densityfrac = mm.getTemplated2DField<double, DataLayout::MAT_DOM, BSet>("Densityfrac");
  //auto Volfrac = mm.getTemplated2DField<double, DataLayout::MAT_DOM, BSet>("Volfrac");

  //auto Densityfrac = mm.get2dFieldAsSlamBivarMap<double,BSet>("Densityfrac");
  //auto Volfrac = mm.get2dFieldAsSlamBivarMap<double,BSet>("Volfrac");

  auto& Vol = mm.get1dField<double>("Vol");

  std::vector<double> Density_average(ncells, 0.0);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    for (auto& v : Density_average)
      v = 0.0;

    timer.start();

    for (int m = 0; m < nmats; ++m)
    {
      auto Densityfrac_row = Densityfrac(m);
      auto Volfrac_row = Volfrac(m);
      for (int j = 0; j < Volfrac_row.size(); ++j)
      {
        Density_average[Densityfrac_row.index(j)] += Densityfrac_row(j) * Volfrac_row(j);
      }
    }
    for (int ic = 0; ic < ncells; ++ic)
    {
      Density_average[ic] /= Vol[ic];
    }

    timer.record();
    data_checker.check(Density_average);
  }
  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::avg_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_submap_templated, act_perf);

  SLIC_INFO("Average Density                      compute time is "
    << act_perf << " secs\n");
}

//    Average density - Material-Dominant
//    MultiMat - Submap
template<typename BSetType>
void average_density_mat_dom_mm_submap(MultiMat& mm) {
  SLIC_INFO("-- Averaging Density mat-dominant using MultiMat Submap -- templated on concrete bivariate set");

  mm.convertLayoutToMaterialDominant();

  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
                                << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::MAT_DOM);

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();

  auto Densityfrac = mm.getTemplated2DField<double, DataLayout::MAT_DOM, BSetType>("Densityfrac");
  auto Volfrac = mm.getTemplated2DField<double, DataLayout::MAT_DOM, BSetType>("Volfrac");
  
  //auto Densityfrac = mm.get2dFieldAsSlamBivarMap<double,BSetType>("Densityfrac");
  //auto Volfrac = mm.get2dFieldAsSlamBivarMap<double,BSetType>("Volfrac");
  
  auto& Vol = mm.get1dField<double>("Vol");

  std::vector<double> Density_average(ncells, 0.0);

  timer.reset();

  for (int iter = 0 ; iter < ITERMAX ; ++iter)
  {
    for (auto& v : Density_average)
      v = 0.0;

    timer.start();

    for (int m = 0 ; m < nmats ; ++m)
    {
      auto Densityfrac_row = Densityfrac(m);
      auto Volfrac_row = Volfrac(m);
      for (int j = 0 ; j < Volfrac_row.size() ; ++j)
      {
        Density_average[Densityfrac_row.index(j)] += Densityfrac_row(j) * Volfrac_row(j);
      }
    }
    for (int ic = 0 ; ic < ncells ; ++ic)
    {
      Density_average[ic] /= Vol[ic];
    }

    timer.record();
    data_checker.check(Density_average);
  }
  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::avg_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_submap_templated, act_perf);

  SLIC_INFO("Average Density                      compute time is "
            << act_perf << " secs\n");
}

//    Average density - Material-Dominant
//    MultiMat - IndexArray
void average_density_mat_dom_mm_idxarray(MultiMat& mm) {
  SLIC_INFO("-- Averaging Density mat-dominant using MultiMat Index Array --");

  mm.convertLayoutToMaterialDominant();

  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::MAT_DOM);

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double> & Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double> & Volfrac = mm.get2dField<double>("Volfrac");
  MultiMat::Field1D<double> & Vol = mm.get1dField<double>("Vol");

  std::vector<double> Density_average(ncells, 0.0);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    for (auto& v : Density_average)
      v = 0.0;

    timer.start();

    for (int m = 0; m < nmats; ++m)
    {
      auto idxSet = mm.getIndexingSetOfMat(m, mm.getFieldSparsityLayout(0));
      auto cellId = mm.getCellContainingMat(m);
      for (int j = 0; j < idxSet.size(); ++j)
      {
        Density_average[cellId[j]] += Densityfrac[idxSet[j]] * Volfrac[idxSet[j]];
      }
    }
    for (int ic = 0; ic < ncells; ++ic)
    {
      Density_average[ic] /= Vol[ic];
    }

    timer.record();
    data_checker.check(Density_average);
  }
  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::avg_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_idxarray, act_perf);
  SLIC_INFO("Average Density                      compute time is "
    << act_perf << " secs\n");
}

//    Average density - Material-Dominant
//    MultiMat - Iterator
void average_density_mat_dom_mm_iter(MultiMat& mm) {
  SLIC_INFO("-- Averaging Density mat-dominant using MultiMat Iterator --");

  mm.convertLayoutToMaterialDominant();

  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::MAT_DOM);

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double> & Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double> & Volfrac = mm.get2dField<double>("Volfrac");
  MultiMat::Field1D<double> & Vol = mm.get1dField<double>("Vol");

  std::vector<double> Density_average(ncells, 0.0);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    for (auto& v : Density_average)
      v = 0.0;

    timer.start();

    for (int m = 0; m < nmats; ++m)
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
    for (int ic = 0; ic < ncells; ++ic)
    {
      Density_average[ic] /= Vol[ic];
    }

    timer.record();
    data_checker.check(Density_average);
  }
  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::avg_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_iter, act_perf);

  SLIC_INFO("Average Density                      compute time is "
    << act_perf << " secs\n");
}


//    Average density - Material-Dominant 
//    MultiMat - Flat Iterator
void average_density_mat_dom_mm_flatiter(MultiMat& mm) {
  SLIC_INFO("-- Averaging Density mat-dominant using MultiMat Flat Iter --");

  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));

  int ncells = mm.getNumberOfCells();
  MultiMat::Field2D<double> & Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double> & Volfrac = mm.get2dField<double>("Volfrac");
  MultiMat::Field1D<double> & Vol = mm.get1dField<double>("Vol");

  std::vector<double> Density_average(ncells, 0.0);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    for (auto& v : Density_average)
      v = 0.0;

    timer.start();

    auto DensityIter = Densityfrac.begin();
    auto VolfracIter = Volfrac.begin();

    auto DensityIterEnd = Densityfrac.end();
    for (; DensityIter != DensityIterEnd; ++DensityIter, ++VolfracIter)
    {
      Density_average[DensityIter.secondIndex()] += *DensityIter * *VolfracIter;
    }

    for (int ic = 0; ic < ncells; ++ic)
    {
      Density_average[ic] /= Vol[ic];
    }

    timer.record();
    data_checker.check(Density_average);
  }
  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::avg_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_flatiter, act_perf);
  SLIC_INFO("Average Density                      compute time is "
    << act_perf << " secs\n");
}


////////////////////// Other things /////////////////////////

//    Average density with if - Cell-Dominant Full Matrix
//      Robey's 
void average_density_cell_dom_with_if(Robey_data& data)
{
  int ncells = data.ncells;
  int nmats = data.nmats;
  std::vector<double>& Volfrac = data.Volfrac;
  std::vector<double>& Densityfrac = data.Densityfrac;
  std::vector<double>& Vol = data.Vol;

  SLIC_INFO("-- Averaging Density with if --");

  std::vector<double> Density_average(ncells, 0.0);
  timer.reset();

  for (int iter = 0 ; iter < ITERMAX ; ++iter)
  {
    timer.start();

    for (int ic = 0 ; ic < ncells ; ++ic)
    {
      double density_ave = 0.0;
      for (int m = 0 ; m < nmats ; ++m)
      {
        if (Volfrac[ic*nmats + m] > 0.0)
        {
          density_ave += Densityfrac[ic*nmats + m] * Volfrac[ic*nmats + m];
        }
      }
      Density_average[ic] = density_ave / Vol[ic];
    }

    timer.record();
  }

  double act_perf = timer.get_median();
  SLIC_INFO("Average Density of frac with if            compute time is "
            << act_perf << " secs\n");

}

//    Average density with if - Cell-Dominant
//      MultiMat 
void average_density_cell_dom_with_if_mm(MultiMat&) {
  //SLIC_INFO("-- Averaging Density with if using MultiMat --");
  //this doesn't really make sense for a MultiMat stored in sparse layout...
}


////////////////////// Calculate pressure - Cell Dominant /////////////////////////

//   Calculate pressure using ideal gas law - Cell-Dominant Full Matrix 
//     Robey's
void calculate_pressure_cell_dom_full(Robey_data& data)
{
  int ncells = data.ncells;
  int nmats = data.nmats;
  std::vector<double>& Volfrac = data.Volfrac;
  std::vector<double>& Densityfrac = data.Densityfrac;
  std::vector<double>& Temperaturefrac = data.Temperaturefrac;
  std::vector<double>& nmatconsts = data.nmatconsts;

  SLIC_INFO("-- Calculating pressure Cell-Dominant Full Matrix array access--");
  std::vector<double> Pressurefrac(ncells*nmats, 0);

  timer.reset();

  for (int iter = 0 ; iter< ITERMAX ; ++iter)
  {
    timer.start();

    for (int ic = 0 ; ic < ncells ; ++ic)
    {
      for (int m = 0 ; m < nmats ; ++m)
      {
        const auto vf = Volfrac[ic*nmats + m];
        
        Pressurefrac[ic*nmats + m] = (vf > 0.)
            ? (nmatconsts[m] * Densityfrac[ic*nmats + m] * Temperaturefrac[ic*nmats + m]) / vf
            : 0.0;
      }
    }

    timer.record();
    data_checker.check(Pressurefrac);
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::pressure_calc, DataLayout::CELL_DOM,
     SparsityLayout::DENSE, Result_Store::method_csr, act_perf);
  SLIC_INFO("Pressure Calculation with if           compute time is " 
      << act_perf << " secs\n");
}

//   Calculate pressure using ideal gas law - Cell-Dominant compact
//     Robey's
void calculate_pressure_cell_dom_compact(Robey_data& data)
{
  int ncells = data.ncells;
  std::vector<double>& Volfrac = data.Volfrac_sparse;
  std::vector<double>& Densityfrac = data.Densityfrac_sparse;
  std::vector<double>& Temperaturefrac = data.Temperaturefrac_sparse;

  std::vector<double>& nmatconsts = data.nmatconsts;
  std::vector<int>& begin_idx = data.begin_idx;
  std::vector<int>& mat_id = data.col_idx;

  SLIC_INFO("-- Calculating pressure Cell-Dominant compact array access--");
  std::vector<double> Pressurefrac(begin_idx.back(), 0);

  timer.reset();

  for (int iter = 0; iter< ITERMAX; ++iter)
  {
    timer.start();

    for (int ic = 0; ic < ncells; ++ic)
    {
      for (int ii = begin_idx[ic]; ii < begin_idx[ic + 1]; ++ii)
      {
        const int m = mat_id[ii];
        const auto vf = Volfrac[ii];

        Pressurefrac[ii] = (vf > 0.)
            ? (nmatconsts[m] * Densityfrac[ii] * Temperaturefrac[ii]) / vf
            :  0.0;
      }
    }

    timer.record();
    data_checker.check(Pressurefrac);
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::pressure_calc, DataLayout::CELL_DOM,
      SparsityLayout::SPARSE, Result_Store::method_csr, act_perf);
  SLIC_INFO("Pressure Calculation with if           compute time is "
      << act_perf << " secs\n");
}


//   Calculate pressure using ideal gas law - Cell-Dominant
//     MultiMat - Direct Access
void calculate_pressure_cell_dom_full_mm_direct(MultiMat& mm)
{
  SLIC_INFO("-- Calculating pressure, using MultiMat Direct Access --");
  mm.convertLayoutToCellDominant();
  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::CELL_DOM);
  SLIC_ASSERT(mm.getFieldSparsityLayout(0) == SparsityLayout::DENSE);

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double> & Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double> & Volfrac = mm.get2dField<double>("Volfrac");
  MultiMat::Field2D<double> & Temperaturefrac =
    mm.get2dField<double>("Tempfrac");
  MultiMat::Field1D<double> & nmatconsts = mm.get1dField<double>("nmatconsts");

  std::vector<double> Pressurefrac(ncells * nmats, 0.0);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    timer.start();

    for (int ic = 0; ic < ncells; ++ic)
    {
      for (int m = 0; m < nmats; ++m)
      {
        if (*Volfrac.findValue(ic,m) > 0.)
        {
          Pressurefrac[ic * nmats + m] = (nmatconsts[m] * *Densityfrac.findValue(ic, m) *
            *Temperaturefrac.findValue(ic, m)) / *Volfrac.findValue(ic, m);
        }
        else {
          Pressurefrac[ic * nmats + m] = 0.0;
        }
      }
    }

    timer.record();
    data_checker.check(Pressurefrac);
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::pressure_calc, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_direct, act_perf);
  SLIC_INFO("Pressure Calculation with if           compute time is "
    << act_perf << " secs\n");
}

//   Calculate pressure using ideal gas law - Cell-Dominant
//     MultiMat - Submap
void calculate_pressure_cell_dom_mm_submap(MultiMat& mm)
{
  SLIC_INFO("-- Calculating pressure, using MultiMat Submap --");
  mm.convertLayoutToCellDominant();
  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
                                << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::CELL_DOM);

  int ncells = mm.getNumberOfCells();
  auto& Densityfrac = mm.get2dField<double>("Densityfrac");
  auto& Volfrac = mm.get2dField<double>("Volfrac");
  auto& Temperaturefrac = mm.get2dField<double>("Tempfrac");
  auto& nmatconsts = mm.get1dField<double>("nmatconsts");
  
  // Store result in a dense field
  /*
  using DenseFieldSet = typename MultiMat::ProductSetType;
  using DenseField = MultiMat::Field2D<double, DenseFieldSet>;
  DenseField Pressurefrac(mm.getDense2dFieldSet());
  */
  auto Pressurefrac = mm.get2dField<double>("Pressurefrac");

  timer.reset();

  for (int iter = 0 ; iter< ITERMAX ; ++iter)
  {
    timer.start();

    for (int ic = 0 ; ic < ncells ; ++ic)
    {
      const auto Densityfrac_row = Densityfrac(ic);
      const auto Volfrac_row = Volfrac(ic);
      const auto Tempfrac_row = Temperaturefrac(ic);
      auto Pressurefrac_row = Pressurefrac(ic);

      const auto cellMats = Volfrac_row.size();
      for (int j = 0; j < cellMats; ++j)
      {
        const int m = Volfrac_row.index(j);
        const auto vf = Volfrac_row(j);

        Pressurefrac_row(j) = vf > 0.
            ? (nmatconsts[m] * Densityfrac_row(j) * Tempfrac_row(j)) / vf
            : 0.;
      }
    }

    timer.record();
    data_checker.check(Pressurefrac.getMap()->data());
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::pressure_calc, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_submap, act_perf);
  SLIC_INFO("Pressure Calculation with if           compute time is "
    << act_perf << " secs\n");
}

template<typename BSet>
void calculate_pressure_cell_dom_mm_submap(MultiMat& mm)
{
  SLIC_INFO("-- Calculating pressure, using MultiMat Submap \t\t-- templated on BSet");
  mm.convertLayoutToCellDominant();
  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
      << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::CELL_DOM);

  int ncells = mm.getNumberOfCells();
  
  //auto Densityfrac = mm.getTemplated2DField<double, DataLayout::CELL_DOM, BSet>("Densityfrac");
  //auto Volfrac = mm.getTemplated2DField<double, DataLayout::CELL_DOM, BSet>("Volfrac");
  //auto Temperaturefrac = mm.getTemplated2DField<double, DataLayout::CELL_DOM, BSet>("Tempfrac");

  auto Densityfrac = mm.get2dFieldAsSlamBivarMap<double,BSet>("Densityfrac");
  auto Volfrac = mm.get2dFieldAsSlamBivarMap<double,BSet>("Volfrac");
  auto Temperaturefrac = mm.get2dFieldAsSlamBivarMap<double,BSet>("Tempfrac");
  
  auto& nmatconsts = mm.get1dField<double>("nmatconsts");

  //auto Pressurefrac = mm.getTemplated2DField<double, DataLayout::CELL_DOM, BSet>("Pressurefrac");
  auto Pressurefrac = mm.get2dFieldAsSlamBivarMap<double, BSet>("Pressurefrac");

  timer.reset();
  
  for (int iter = 0; iter< ITERMAX; ++iter)
  {
    timer.start();
    
    for (int ic = 0; ic < ncells; ++ic)
    {
      const auto Densityfrac_row = Densityfrac(ic);
      const auto Volfrac_row = Volfrac(ic);
      const auto Tempfrac_row = Temperaturefrac(ic);
      auto Pressurefrac_row = Pressurefrac(ic);

      const auto cellMats = Volfrac_row.size();
      for (int j = 0; j < cellMats; ++j)
      {
        const int m = Volfrac_row.index(j);
        const auto vf = Volfrac_row(j);

        //Pressurefrac_row(m) 
        Pressurefrac_row(j) = (vf > 0.)
            ? (nmatconsts[m] * Densityfrac_row(j) * Tempfrac_row(j)) / vf
            : 0.;
      }
    }

    timer.record();
    data_checker.check(Pressurefrac.getMap()->data());
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::pressure_calc, mm.getFieldDataLayout(0),
      mm.getFieldSparsityLayout(0), Result_Store::mm_submap_templated, act_perf);
  SLIC_INFO("Pressure Calculation with if           compute time is "
      << act_perf << " secs\n");
}

template<typename BSet, typename FieldTypeT>
void calculate_pressure_cell_dom_mm_submap(MultiMat& mm)
{
  SLIC_INFO("-- Calculating pressure, using MultiMat Submap \t\t-- templated on BSet");
  mm.convertLayoutToCellDominant();
  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::CELL_DOM);

  int ncells = mm.getNumberOfCells();

  //pick which get field function to use
  GetFieldFunType<FieldTypeT> getFieldFun = getfieldfun<FieldTypeT, BSet>();

  FieldTypeT Densityfrac = (mm.*getFieldFun)("Densityfrac");
  FieldTypeT Volfrac = (mm.*getFieldFun)("Volfrac");
  FieldTypeT Temperaturefrac = (mm.*getFieldFun)("Tempfrac");
  FieldTypeT Pressurefrac = (mm.*getFieldFun)("Pressurefrac");

  
  auto& nmatconsts = mm.get1dField<double>("nmatconsts");

  //auto Pressurefrac = mm.getTemplated2DField<double, DataLayout::CELL_DOM, BSet>("Pressurefrac");
  //auto Pressurefrac = mm.get2dFieldAsSlamBivarMap<double, BSet>("Pressurefrac");

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    timer.start();

    for (int ic = 0; ic < ncells; ++ic)
    {
      const auto Densityfrac_row = Densityfrac(ic);
      const auto Volfrac_row = Volfrac(ic);
      const auto Tempfrac_row = Temperaturefrac(ic);
      auto Pressurefrac_row = Pressurefrac(ic);

      const auto cellMats = Volfrac_row.size();
      for (int j = 0; j < cellMats; ++j)
      {
        const int m = Volfrac_row.index(j);
        const auto vf = Volfrac_row(j);

        //Pressurefrac_row(m) 
        Pressurefrac_row(j) = (vf > 0.)
          ? (nmatconsts[m] * Densityfrac_row(j) * Tempfrac_row(j)) / vf
          : 0.;
      }
    }

    timer.record();
    data_checker.check(Pressurefrac.getMap()->data());
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::pressure_calc, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_submap_templated, act_perf);
  SLIC_INFO("Pressure Calculation with if           compute time is "
    << act_perf << " secs\n");
}

//   Calculate pressure using ideal gas law - Cell-Dominant
//     MultiMat - Iterator
void calculate_pressure_cell_dom_full_mm_iter(MultiMat& mm)
{
  SLIC_INFO("-- Calculating pressure, using MultiMat Iterator --");
  mm.convertLayoutToCellDominant();
  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::CELL_DOM);
  SLIC_ASSERT(mm.getFieldSparsityLayout(0) == SparsityLayout::DENSE);

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double> & Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double> & Volfrac = mm.get2dField<double>("Volfrac");
  MultiMat::Field2D<double> & Temperaturefrac =
    mm.get2dField<double>("Tempfrac");
  MultiMat::Field1D<double> & nmatconsts = mm.get1dField<double>("nmatconsts");

  std::vector<double> Pressurefrac(ncells * nmats, 0.0);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    timer.start();

    for (int ic = 0; ic < ncells; ++ic)
    {
      auto DensityfracIter = Densityfrac.begin(ic);
      auto VolfracIter = Volfrac.begin(ic);
      auto TempfracIter = Temperaturefrac.begin(ic);

      auto VolfracIterEnd = Volfrac.end(ic);
      for (; VolfracIter != VolfracIterEnd; ++VolfracIter, ++DensityfracIter, ++TempfracIter)
      {
        int m = VolfracIter.index();
        if (*VolfracIter > 0.)
        {
          Pressurefrac[ic * nmats + m] = (nmatconsts[m] * *DensityfracIter *
            *TempfracIter) / *VolfracIter;
        }
        else {
          Pressurefrac[ic * nmats + m] = 0;
        }
      }
    }

    timer.record();
    data_checker.check(Pressurefrac);
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::pressure_calc, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_iter, act_perf);
  SLIC_INFO("Pressure Calculation with if           compute time is "
    << act_perf << " secs\n");
}

//   Calculate pressure using ideal gas law - Cell-Dominant
//     MultiMat - Submap
void calculate_pressure_cell_dom_full_mm_flatiter(MultiMat& mm)
{
  SLIC_INFO("-- Calculating pressure, using MultiMat Flat Iterator --");
  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldSparsityLayout(0) == SparsityLayout::DENSE);
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::CELL_DOM);

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double> & Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double> & Volfrac = mm.get2dField<double>("Volfrac");
  MultiMat::Field2D<double> & Temperaturefrac =
    mm.get2dField<double>("Tempfrac");
  MultiMat::Field1D<double> & nmatconsts = mm.get1dField<double>("nmatconsts");

  std::vector<double> Pressurefrac(ncells * nmats, 0.0);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    timer.start();

    auto DensityfracIter = Densityfrac.begin();
    auto VolfracIter = Volfrac.begin();
    auto TempfracIter = Temperaturefrac.begin();

    int ii = 0;
    auto VolfracIterEnd = Volfrac.end();
    for (; VolfracIter != VolfracIterEnd; ++VolfracIter, ++DensityfracIter, ++TempfracIter, ++ii)
    {
      int m = VolfracIter.secondIndex();
      if (*VolfracIter > 0.)
      {
        Pressurefrac[ii] = (nmatconsts[m] * *DensityfracIter *
          *TempfracIter) / *VolfracIter;
      }
      else {
        Pressurefrac[ii] = 0;
      }
    }

    timer.record();
    data_checker.check(Pressurefrac);
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::pressure_calc, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_flatiter, act_perf);
  SLIC_INFO("Pressure Calculation with if           compute time is "
    << act_perf << " secs\n");
}

////////////////////// Calculate pressure - Material Dominant /////////////////////////


//   Calculate pressure using ideal gas law - Material-Dominant Full Matrix 
//     Robey's
void calculate_pressure_mat_dom_full(Robey_data& data)
{
  int ncells = data.ncells;
  int nmats = data.nmats;
  std::vector<double>& Volfrac = data.Volfrac;
  std::vector<double>& Densityfrac = data.Densityfrac;
  std::vector<double>& Temperaturefrac = data.Temperaturefrac;
  std::vector<double>& nmatconsts = data.nmatconsts;

  SLIC_INFO("-- Calculating pressure Material-Dominant Full Matrix array access--");
  std::vector<double> Pressurefrac(ncells*nmats, 0);

  timer.reset();

  for (int iter = 0 ; iter< ITERMAX ; ++iter)
  {
    timer.start();
    for (int m = 0; m < nmats; ++m) {
      for (int ic = 0; ic < ncells; ++ic) {
        if (Volfrac[m*ncells+ic] > 0.0) {
          Pressurefrac[m * ncells + ic] = (nmatconsts[m] * Densityfrac[m * ncells + ic] *
            Temperaturefrac[m * ncells + ic]) / Volfrac[m * ncells + ic];
        }
        else {
          Pressurefrac[m * ncells + ic] = 0.0;
        }
      }
    }

    timer.record();
    data_checker.check(Pressurefrac);
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::pressure_calc, DataLayout::MAT_DOM,
     SparsityLayout::DENSE, Result_Store::method_csr, act_perf);
  SLIC_INFO("Pressure Calculation with if           compute time is " 
      << act_perf << " secs\n");
}

//   Calculate pressure using ideal gas law - Material-Dominant compact Matrix 
//     Robey's
void calculate_pressure_mat_dom_compact(Robey_data& data)
{
  int ncells = data.ncells;
  int nmats = data.nmats;
  std::vector<double>& Volfrac = data.Volfrac_sparse;
  std::vector<double>& Densityfrac = data.Densityfrac_sparse;
  std::vector<double>& Temperaturefrac = data.Temperaturefrac_sparse;
  std::vector<double>& nmatconsts = data.nmatconsts;
  std::vector<int>& begin_idx = data.begin_idx;
  std::vector<int>& cell_id = data.col_idx;

  SLIC_INFO("-- Calculating pressure Material-Dominant compact array access--");
  std::vector<double> Pressurefrac(ncells*nmats, 0);

  timer.reset();

  for (int iter = 0; iter< ITERMAX; ++iter)
  {
    timer.start();
    for (int m = 0; m < nmats; ++m) 
    {
      const auto matconst = nmatconsts[m];

      for (int ii = begin_idx[m]; ii < begin_idx[m + 1]; ++ii)
      { 
        const int ci = cell_id[ii];
        const auto vf = Volfrac[ii];

        Pressurefrac[m * ncells + ci] = (vf > 0.0)
          ? ( matconst * Densityfrac[ii] * Temperaturefrac[ii]) / vf
          : 0.;
      }
    }

    timer.record();
    data_checker.check(Pressurefrac);
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::pressure_calc, DataLayout::MAT_DOM,
        SparsityLayout::SPARSE, Result_Store::method_csr, act_perf);
  SLIC_INFO("Pressure Calculation with if           compute time is "
        << act_perf << " secs\n");
}


//   Calculate pressure using ideal gas law - Material-Dominant
//     MultiMat - Direct Access
void calculate_pressure_mat_dom_full_mm_direct(MultiMat& mm)
{
  SLIC_INFO("-- Calculating pressure, using MultiMat Direct Access --");
  mm.convertLayoutToMaterialDominant();
  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::MAT_DOM);
  SLIC_ASSERT(mm.getFieldSparsityLayout(0) == SparsityLayout::DENSE);

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double> & Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double> & Volfrac = mm.get2dField<double>("Volfrac");
  MultiMat::Field2D<double> & Temperaturefrac =
    mm.get2dField<double>("Tempfrac");
  MultiMat::Field1D<double> & nmatconsts = mm.get1dField<double>("nmatconsts");

  std::vector<double> Pressurefrac(ncells * nmats, 0.0);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    timer.start();

    for (int m = 0; m < nmats; ++m) {
      for (int ic = 0; ic < ncells; ++ic) {
        if (*Volfrac.findValue(m, ic) > 0.)
        {
          Pressurefrac[m*ncells+ic] = (nmatconsts[m] * *Densityfrac.findValue(m,ic) *
            *Temperaturefrac.findValue(m, ic)) / *Volfrac.findValue(m, ic);
        }
        else {
          Pressurefrac[m * ncells + ic] = 0.0;
        }
      }
    }

    timer.record();
    data_checker.check(Pressurefrac);
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::pressure_calc, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_direct, act_perf);
  SLIC_INFO("Pressure Calculation with if           compute time is "
    << act_perf << " secs\n");
}

//   Calculate pressure using ideal gas law - Material-Dominant
//     MultiMat - Submap
void calculate_pressure_mat_dom_mm_submap(MultiMat& mm)
{
  SLIC_INFO("-- Calculating pressure, using MultiMat Submap --");
  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
                                << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::MAT_DOM);

  int nmats = mm.getNumberOfMaterials();
  auto& Densityfrac = mm.get2dField<double>("Densityfrac");
  auto& Volfrac = mm.get2dField<double>("Volfrac");
  auto& Temperaturefrac = mm.get2dField<double>("Tempfrac");
  auto& nmatconsts = mm.get1dField<double>("nmatconsts");
  
  // Store result in a dense field
  using DenseFieldSet = typename MultiMat::ProductSetType;
  using DenseField = MultiMat::Field2D<double, DenseFieldSet>;

  DenseField Pressurefrac(mm, mm.getDense2dFieldSet(DataLayout::MAT_DOM));

  timer.reset();

  for (int iter = 0 ; iter< ITERMAX ; ++iter)
  {
    timer.start();

    for (int m = 0; m < nmats; ++m) 
    {
      const auto Densityfrac_row = Densityfrac(m);
      const auto Volfrac_row = Volfrac(m);
      const auto Tempfrac_row = Temperaturefrac(m);
      auto Pressurefrac_row = Pressurefrac(m);

      const auto mprop = nmatconsts[m];

      const int matCells = Volfrac_row.size();
      for (int j = 0; j < matCells; ++j)
      {
        const int ic = Volfrac_row.index(j);
        const auto vf = Volfrac_row(j);

        Pressurefrac_row(ic) = (vf > 0.)
            ? (mprop * Densityfrac_row(j) * Tempfrac_row(j)) / vf
            : 0.;
      }
    }

    timer.record();
    data_checker.check(Pressurefrac.getMap()->data());
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::pressure_calc, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_submap, act_perf);
  SLIC_INFO("Pressure Calculation with if           compute time is "
    << act_perf << " secs\n");
}

template<typename BSet>
void calculate_pressure_mat_dom_mm_submap(MultiMat& mm)
{
  SLIC_INFO("-- Calculating pressure, using MultiMat Submap \t\t-- using templated BSet");
  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
      << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::MAT_DOM);

  int nmats = mm.getNumberOfMaterials();
  
  //auto Densityfrac = mm.getTemplated2DField<double, DataLayout::MAT_DOM, BSet>("Densityfrac");
  //auto Volfrac = mm.getTemplated2DField<double, DataLayout::MAT_DOM, BSet>("Volfrac");
  //auto Temperaturefrac = mm.getTemplated2DField<double, DataLayout::MAT_DOM, BSet>("Tempfrac");

  auto Densityfrac = mm.get2dFieldAsSlamBivarMap<double,BSet>("Densityfrac");
  auto Volfrac = mm.get2dFieldAsSlamBivarMap<double,BSet>("Volfrac");
  auto Temperaturefrac = mm.get2dField<double>("Tempfrac");

  MultiMat::Field1D<double>& nmatconsts = mm.get1dField<double>("nmatconsts");

  // Store result in a dense field
  using DenseFieldSet = typename MultiMat::ProductSetType;
  using DenseField = MultiMat::Field2D<double, DenseFieldSet>;
  DenseField Pressurefrac(mm, mm.getDense2dFieldSet(DataLayout::MAT_DOM));

  //this field needs to be dense, but the other fields can be sparse, so not this:
  //auto Pressurefrac = mm.getTemplated2DField<double, DataLayout::MAT_DOM, BSet>("Pressurefrac");

  timer.reset();

  for (int iter = 0; iter< ITERMAX; ++iter)
  {
    timer.start();

    for (int m = 0; m < nmats; ++m)
    {
      const auto Densityfrac_row = Densityfrac(m);
      const auto Volfrac_row = Volfrac(m);
      const auto Tempfrac_row = Temperaturefrac(m);
      auto Pressurefrac_row = Pressurefrac(m);

      const auto mprop = nmatconsts[m];

      const int matCells = Volfrac_row.size();
      for (int j = 0; j < matCells; ++j)
      {
        const int ic = Volfrac_row.index(j);
        const auto vf = Volfrac_row(j);

        Pressurefrac_row(ic) = (vf > 0)
            ? (mprop * Densityfrac_row(j) * Tempfrac_row(j)) / vf
            : 0;
      }
    }

    timer.record();
    data_checker.check(Pressurefrac.getMap()->data());
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::pressure_calc, mm.getFieldDataLayout(0),
      mm.getFieldSparsityLayout(0), Result_Store::mm_submap_templated, act_perf);
  SLIC_INFO("Pressure Calculation with if           compute time is "
      << act_perf << " secs\n");
}


template<typename BSet, typename FieldTypeT>
void calculate_pressure_mat_dom_mm_submap(MultiMat& mm)
{
  SLIC_INFO("-- Calculating pressure, using MultiMat Submap \t\t-- using templated BSet");
  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::MAT_DOM);

  int nmats = mm.getNumberOfMaterials();

  //pick which get field function to use
  GetFieldFunType<FieldTypeT> getFieldFun = getfieldfun<FieldTypeT, BSet>();

  FieldTypeT Densityfrac = (mm.*getFieldFun)("Densityfrac");
  FieldTypeT Volfrac = (mm.*getFieldFun)("Volfrac");
  FieldTypeT Temperaturefrac = (mm.*getFieldFun)("Tempfrac");

  MultiMat::Field1D<double>& nmatconsts = mm.get1dField<double>("nmatconsts");

  // Store result in a dense field
  using DenseFieldSet = typename MultiMat::ProductSetType;
  using DenseField = MultiMat::Field2D<double, DenseFieldSet>;
  DenseField Pressurefrac(mm, mm.getDense2dFieldSet(DataLayout::MAT_DOM));

  //this field needs to be dense, but the other fields can be sparse, so not this:
  //auto Pressurefrac = mm.getTemplated2DField<double, DataLayout::MAT_DOM, BSet>("Pressurefrac");

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    timer.start();

    for (int m = 0; m < nmats; ++m)
    {
      const auto Densityfrac_row = Densityfrac(m);
      const auto Volfrac_row = Volfrac(m);
      const auto Tempfrac_row = Temperaturefrac(m);
      auto Pressurefrac_row = Pressurefrac(m);

      const auto mprop = nmatconsts[m];

      const int matCells = Volfrac_row.size();
      for (int j = 0; j < matCells; ++j)
      {
        const int ic = Volfrac_row.index(j);
        const auto vf = Volfrac_row(j);

        Pressurefrac_row(ic) = (vf > 0)
          ? (mprop * Densityfrac_row(j) * Tempfrac_row(j)) / vf
          : 0;
      }
    }

    timer.record();
    data_checker.check(Pressurefrac.getMap()->data());
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::pressure_calc, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_submap_templated, act_perf);
  SLIC_INFO("Pressure Calculation with if           compute time is "
    << act_perf << " secs\n");
}


//   Calculate pressure using ideal gas law - Material-Dominant
//     MultiMat - Iterator
void calculate_pressure_mat_dom_full_mm_iter(MultiMat& mm)
{
  SLIC_INFO("-- Calculating pressure, using MultiMat Iter --");
  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::MAT_DOM);
  SLIC_ASSERT(mm.getFieldSparsityLayout(0) == SparsityLayout::DENSE);

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double> & Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double> & Volfrac = mm.get2dField<double>("Volfrac");
  MultiMat::Field2D<double> & Temperaturefrac =
    mm.get2dField<double>("Tempfrac");
  MultiMat::Field1D<double> & nmatconsts = mm.get1dField<double>("nmatconsts");

  std::vector<double> Pressurefrac(ncells * nmats, 0.0);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    timer.start();

    for (int m = 0; m < nmats; ++m)
    {
      auto DensityfracIter = Densityfrac.begin(m);
      auto VolfracIter = Volfrac.begin(m);
      auto TempfracIter = Temperaturefrac.begin(m);

      auto VolfracIterEnd = Volfrac.end(m);
      for (; VolfracIter != VolfracIterEnd; ++VolfracIter, ++DensityfracIter, ++TempfracIter)
      {
        int ic = VolfracIter.index();
        if (*VolfracIter > 0.)
        {
          Pressurefrac[m * ncells + ic] = (nmatconsts[m] * *DensityfracIter *
            *TempfracIter) / *VolfracIter;
        }
        else {
          Pressurefrac[m * ncells + ic] = 0;
        }
      }
    }

    timer.record();
    data_checker.check(Pressurefrac);
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::pressure_calc, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_iter, act_perf);
  SLIC_INFO("Pressure Calculation with if           compute time is "
    << act_perf << " secs\n");
}

//   Calculate pressure using ideal gas law - Cell-Dominant
//     MultiMat - Submap
void calculate_pressure_mat_dom_full_mm_flatiter(MultiMat& mm)
{
  SLIC_INFO("-- Calculating pressure, using MultiMat Flat Iterator --");
  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldSparsityLayout(0) == SparsityLayout::DENSE);
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::MAT_DOM);

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double> & Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double> & Volfrac = mm.get2dField<double>("Volfrac");
  MultiMat::Field2D<double> & Temperaturefrac =
    mm.get2dField<double>("Tempfrac");
  MultiMat::Field1D<double> & nmatconsts = mm.get1dField<double>("nmatconsts");

  std::vector<double> Pressurefrac(ncells * nmats, 0.0);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    timer.start();

    auto DensityfracIter = Densityfrac.begin();
    auto VolfracIter = Volfrac.begin();
    auto TempfracIter = Temperaturefrac.begin();

    int ii = 0;
    auto VolfracIterEnd = Volfrac.end();
    for (; VolfracIter != VolfracIterEnd; ++VolfracIter, ++DensityfracIter, ++TempfracIter, ++ii)
    {
      int m = VolfracIter.firstIndex();
      if (*VolfracIter > 0.)
      {
        Pressurefrac[ii] = (nmatconsts[m] * *DensityfracIter *
          *TempfracIter) / *VolfracIter;
      }
      else {
        Pressurefrac[ii] = 0;
      }
    }

    timer.record();
    data_checker.check(Pressurefrac);
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::pressure_calc, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_flatiter, act_perf);
  SLIC_INFO("Pressure Calculation with if           compute time is "
    << act_perf << " secs\n");
}


////////////////////// Avg Neighbor - Cell Dominant /////////////////////////

//    Average density over neighbor cells
//      Robey's - Cell-Dominant Full Matrix
void average_density_over_nbr_cell_dom_full(Robey_data& data)
{
  int ncells = data.ncells;
  int nmats = data.nmats;
  std::vector<double>& Volfrac = data.Volfrac;
  std::vector<double>& Densityfrac = data.Densityfrac;
  const std::vector<double>& cen = data.cen;
  const std::vector<int>& nnbrs = data.nnbrs;
  const std::vector<int>& nbrs = data.nbrs;


  SLIC_INFO("-- Average Density over Neighbors, Cell-Dominant Full Matrix Array Access --");

  std::vector<double> MatDensity_average(ncells * nmats, 0.0);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    for (auto& v : MatDensity_average)
      v = 0.0;

    timer.start();

    for (int ic = 0; ic < ncells; ++ic)
    {
      double xc[2] = { cen[ic*2], cen[ic*2 + 1] };
      int nn = nnbrs[ic];
      const int* cnbrs = &(nbrs[ic * 8]);
      
      for (int m = 0; m < nmats; ++m)
      {
        if (Volfrac[ic * nmats + m] > 0.0)
        {
          int nnm = 0;         // number of nbrs with this material
          double den = 0.0;
          for (int n = 0; n < nn; ++n)
          {
            int jc = cnbrs[n];
            if (Volfrac[jc * nmats + m] > 0.0)
            {
              double dx = xc[0] - cen[jc * 2 + 0];
              double dy = xc[1] - cen[jc * 2 + 1];
              double dsqr = dx * dx + dy * dy;

              den += Densityfrac[jc * nmats + m] / dsqr;
              ++nnm;
            }
          }
          if (nnm > 0)
            MatDensity_average[ic * nmats + m] = den / nnm;
          else
            SLIC_ASSERT(MatDensity_average[ic * nmats + m] == 0.0);
        }
        else
        {
          MatDensity_average[ic * nmats + m] = 0.0;
        }
      }
    }

    timer.record();
    data_checker.check(MatDensity_average);
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::neighbor_density, DataLayout::CELL_DOM,
    SparsityLayout::DENSE, Result_Store::method_csr, act_perf);
  SLIC_INFO("Average Material Density            compute time is "
    << act_perf << " secs\n");

}

//    Average density over neighbor cells
//      Robey's - Cell-Dominant Compact
void average_density_over_nbr_cell_dom_compact(Robey_data& data)
{
  int ncells = data.ncells;
  std::vector<double>& Densityfrac = data.Densityfrac_sparse;
  const std::vector<double>& cen = data.cen;
  const std::vector<int>& nnbrs = data.nnbrs;
  const std::vector<int>& nbrs = data.nbrs;
  std::vector<int>& begin_idx = data.begin_idx;
  std::vector<int>& mat_id = data.col_idx;

  SLIC_INFO("-- Average Density over Neighbors, Cell-Dominant Compact Matrix Array Access --");

  std::vector<double> MatDensity_average(begin_idx.back(), 0.0);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    for (auto& v : MatDensity_average)
      v = 0.0;

    timer.start();

    for (int ic = 0; ic < ncells; ++ic)
    {
      double xc[2] = { cen[ic * 2], cen[ic * 2 + 1] };
      int nn = nnbrs[ic];
      const int* cnbrs = &(nbrs[ic * 8]);

      for (int ii = begin_idx[ic]; ii < begin_idx[ic + 1]; ++ii)
      {
        int m = mat_id[ii];
        double den = 0.0;
        int nnm = 0;         // number of nbrs with this material
        for (int n = 0; n < nn; ++n)
        {
          int jc = cnbrs[n];
          for (int jj = begin_idx[jc]; jj < begin_idx[jc + 1]; ++jj) //loop through all map in cell jc to find m
          {
            if (mat_id[jj] == m)
            {
              double dx = xc[0] - cen[jc * 2 + 0];
              double dy = xc[1] - cen[jc * 2 + 1];
              double dsqr = dx * dx + dy * dy;

              den += Densityfrac[jj] / dsqr;
              ++nnm;
              break;
            }
          }
        }

        if(nnm > 0)
          MatDensity_average[ii] = den / nnm;
        else
          SLIC_ASSERT(MatDensity_average[ii] == 0.0);
      }
    }
    
    timer.record();
    data_checker.check(MatDensity_average);
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::neighbor_density, DataLayout::CELL_DOM,
    SparsityLayout::SPARSE, Result_Store::method_csr, act_perf);
  SLIC_INFO("Average Material Density            compute time is "
    << act_perf << " secs\n");

}



//    Average material density over neighborhood of each cell
//      MultiMat - Dense - Direct Access
void average_density_over_nbr_cell_dom_full_mm_direct(MultiMat& mm, Robey_data& data)
{
  SLIC_INFO("-- Average Density over Neighbors, Cell-Dominant Full Matrix,"
    << " Multimat Direct Access --");
  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::CELL_DOM);
  SLIC_ASSERT(mm.getFieldSparsityLayout(0) == SparsityLayout::DENSE);

  const std::vector<double> & cen = data.cen;
  const std::vector<int> & nnbrs = data.nnbrs;
  const std::vector<int> & nbrs = data.nbrs;

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double> & Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double> & Volfrac = mm.get2dField<double>("Volfrac");

  std::vector<double> MatDensity_average(ncells * nmats, 0);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    for (auto& v : MatDensity_average)
      v = 0.0;

    timer.start();

    for (int ic = 0; ic < ncells; ++ic)
    {
      double xc[2];
      xc[0] = cen[ic * 2]; xc[1] = cen[ic * 2 + 1];
      int nn = nnbrs[ic];
      int cnbrs[8];
      double dsqr[8];

      for (int n = 0; n < nn; ++n)
        cnbrs[n] = nbrs[ic * 8 + n];

      for (int n = 0; n < nn; ++n)
      {
        dsqr[n] = 0.0;
        for (int d = 0; d < 2; ++d) //original condition was  d < 1 ???
        {
          double ddist = (xc[d] - cen[cnbrs[n] * 2 + d]);
          dsqr[n] += ddist * ddist;
        }
      }

      for (int m = 0; m < nmats; ++m)
      {
        if (*Volfrac.findValue(ic,m) > 0.0)  //this check is not needed in sparse layout.
        {
          int nnm = 0;         // number of nbrs with this material
          for (int n = 0; n < nn; ++n)
          {
            int jc = cnbrs[n];
            if (*Volfrac.findValue(jc, m) > 0.0)
            {
              MatDensity_average[ic * nmats + m] += *Densityfrac.findValue(jc,m) / dsqr[n];
              ++nnm;
            }
          }
          if (nnm > 0)
            MatDensity_average[ic * nmats + m] /= nnm;
          else
            SLIC_ASSERT(MatDensity_average[ic * nmats + m] == 0.0);
        }
        else
        {
          MatDensity_average[ic * nmats + m] = 0.0;
        }
      }
    }

    timer.record();
    data_checker.check(MatDensity_average);
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::neighbor_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_direct, act_perf);

  SLIC_INFO("Average Material Density            compute time is "
    << act_perf << " secs\n");

}


//    Average material density over neighborhood of each cell
//      MultiMat - Dense - Submap
void average_density_over_nbr_cell_dom_full_mm_submap(MultiMat& mm, Robey_data& data)
{
  SLIC_INFO("-- Average Density over Neighbors, Cell-Dominant Full Matrix,"
            << " Multimat Submap --");
  mm.convertLayoutToCellDominant();
  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
                                << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::CELL_DOM);
  SLIC_ASSERT(mm.getFieldSparsityLayout(0) == SparsityLayout::DENSE);

  const std::vector<double>& cen = data.cen;
  const std::vector<int>& nnbrs = data.nnbrs;
  const std::vector<int>& nbrs = data.nbrs;

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double>& Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double>& Volfrac = mm.get2dField<double>("Volfrac");

  std::vector<double> MatDensity_average(ncells*nmats, 0);
  
  timer.reset();
  
  for (int iter = 0 ; iter< ITERMAX ; ++iter)
  {
    for (auto& v : MatDensity_average)
      v = 0.0;

    timer.start();

    for (int ic = 0 ; ic < ncells ; ++ic)
    {
      auto Densityfrac_row = Densityfrac(ic);
      auto Volfrac_row = Volfrac(ic);

      double xc[2];
      xc[0] = cen[ic * 2]; xc[1] = cen[ic * 2 + 1];
      int nn = nnbrs[ic];
      int cnbrs[8];
      double dsqr[8];

      for (int n = 0 ; n < nn ; ++n)
        cnbrs[n] = nbrs[ic * 8 + n];

      for (int n = 0 ; n < nn ; ++n)
      {
        dsqr[n] = 0.0;
        for (int d = 0 ; d < 2 ; ++d) //original condition was  d < 1 ???
        {
          double ddist = (xc[d] - cen[cnbrs[n] * 2 + d]);
          dsqr[n] += ddist * ddist;
        }
      }

      for (int m = 0; m < nmats; ++m)
      {
        if (Volfrac_row(m)> 0.0)  //this check is not needed in sparse layout.
        {
          int nnm = 0;         // number of nbrs with this material
          for (int n = 0 ; n < nn ; ++n)
          {
            int jc = cnbrs[n];
            auto Volfrac_jrow = Volfrac(jc);
            if( Volfrac_jrow(m) > 0.0)
            {
              auto Densityfrac_jrow = Densityfrac(jc);
              MatDensity_average[ic*nmats + m] += Densityfrac_jrow(m) / dsqr[n];
              ++nnm;
            }
          }
          if (nnm > 0)
            MatDensity_average[ic*nmats + m] /= nnm;
          else
            SLIC_ASSERT(MatDensity_average[ic * nmats + m] == 0.0);
        }
        else
        {
          MatDensity_average[ic*nmats + m] = 0.0;
        }
      }
    }

    timer.record();
    data_checker.check(MatDensity_average);
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::neighbor_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_submap, act_perf);

  SLIC_INFO("Average Material Density            compute time is "
            << act_perf << " secs\n");

}


//    Average material density over neighborhood of each cell
//      MultiMat - Dense - Submap
template<typename BSet, typename FieldTypeT>
void average_density_over_nbr_cell_dom_full_mm_submap(MultiMat& mm, Robey_data& data)
{
  SLIC_INFO("-- Average Density over Neighbors, Cell-Dominant Full Matrix,"
    << " Multimat Submap -- using templated BSet ");
  mm.convertLayoutToCellDominant();
  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::CELL_DOM);
  SLIC_ASSERT(mm.getFieldSparsityLayout(0) == SparsityLayout::DENSE);

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();

  //pick which get field function to use
  GetFieldFunType<FieldTypeT> getFieldFun = getfieldfun<FieldTypeT, BSet>();

  FieldTypeT Densityfrac = (mm.*getFieldFun)("Densityfrac");
  FieldTypeT Volfrac = (mm.*getFieldFun)("Volfrac");

  //auto Densityfrac = mm.getTemplated2DField<double, DataLayout::CELL_DOM, BSet>("Densityfrac");
  //auto Volfrac = mm.getTemplated2DField<double, DataLayout::CELL_DOM, BSet>("Volfrac");

  //auto Densityfrac = mm.get2dFieldAsSlamBivarMap<double,BSet>("Densityfrac");
  //auto Volfrac = mm.get2dFieldAsSlamBivarMap<double,BSet>("Volfrac");

  // Get the slam relations and maps for the mesh data
  const auto& neighbors = data.slam_neighbors;
  const auto& cen = data.slam_centroids;

  FieldTypeT MatDensity_average = (mm.*getFieldFun)("MatDensityAverage");
  //auto MatDensity_average = mm.get2dFieldAsSlamBivarMap<double, BSet>("MatDensityAverage");
  //auto MatDensity_average = mm.getTemplated2DField<double, DataLayout::CELL_DOM, BSet>("MatDensityAverage");

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    MatDensity_average.clear();

    timer.start();

    for (int ic = 0; ic < ncells; ++ic)
    {
      // Get the set of neighbors
      auto nbrs = neighbors[ic];

      // compute the squared distances to neighbors
      double xc[2] = { cen(ic,0), cen(ic,1) };

      auto Volfrac_row = Volfrac(ic);
      auto AvgMatDensity_row = MatDensity_average(ic);

      for (int m = 0; m < nmats; ++m)
      {
        if (Volfrac_row(m) > 0.0)  //this check is not needed in sparse layout.
        {
          double den = 0.;
          int nnm = 0;             // number of nbrs with this material
          for (int n : nbrs.positions())
          {
            int jc = nbrs[n];
            if (Volfrac(jc, m) > 0.0)
            {
              double dx = xc[0] - cen(jc, 0);
              double dy = xc[1] - cen(jc, 1);
              double dsqr = dx * dx + dy * dy;

              den += Densityfrac(jc, m) / dsqr;
              ++nnm;
            }
          }
          AvgMatDensity_row(m) = nnm > 0 ? den / nnm : 0.;
        }
        else
        {
          AvgMatDensity_row(m) = 0.;
        }
      }
    }

    timer.record();
    data_checker.check(MatDensity_average.getMap()->data());
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::neighbor_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_submap_templated, act_perf);

  SLIC_INFO("Average Material Density            compute time is "
    << act_perf << " secs\n");

}


//    Average material density over neighborhood of each cell
//      MultiMat - Dense - Submap
template<typename BSet>
void average_density_over_nbr_cell_dom_full_mm_submap(MultiMat& mm, Robey_data& data)
{
  SLIC_INFO("-- Average Density over Neighbors, Cell-Dominant Full Matrix,"
            << " Multimat Submap -- using templated BSet ");
  mm.convertLayoutToCellDominant();
  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
                                << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::CELL_DOM);
  SLIC_ASSERT(mm.getFieldSparsityLayout(0) == SparsityLayout::DENSE);

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  
  auto Densityfrac = mm.getTemplated2DField<double, DataLayout::CELL_DOM, BSet>("Densityfrac");
  auto Volfrac = mm.getTemplated2DField<double, DataLayout::CELL_DOM, BSet>("Volfrac");
  
  //auto Densityfrac = mm.get2dFieldAsSlamBivarMap<double,BSet>("Densityfrac");
  //auto Volfrac = mm.get2dFieldAsSlamBivarMap<double,BSet>("Volfrac");

  // Get the slam relations and maps for the mesh data
  const auto& neighbors = data.slam_neighbors;
  const auto& cen = data.slam_centroids;

  auto MatDensity_average = mm.getTemplated2DField<double, DataLayout::CELL_DOM, BSet>("MatDensityAverage");

  timer.reset();

  for (int iter = 0 ; iter< ITERMAX ; ++iter)
  {
    MatDensity_average.clear();

    timer.start();

    for (int ic = 0 ; ic < ncells ; ++ic)
    {
      // Get the set of neighbors
      auto nbrs = neighbors[ic];

      // compute the squared distances to neighbors
      double xc[2] = {cen(ic,0), cen(ic,1)};

      auto Volfrac_row = Volfrac(ic);
      auto AvgMatDensity_row = MatDensity_average(ic);

      for (int m = 0; m < nmats; ++m)
      {
        if (Volfrac_row(m) > 0.0)  //this check is not needed in sparse layout.
        {
          double den = 0.;
          int nnm = 0;             // number of nbrs with this material
          for (int n: nbrs.positions())
          {
            int jc = nbrs[n];
            if( Volfrac(jc, m) > 0.0)
            {
              double dx = xc[0] - cen(jc, 0);
              double dy = xc[1] - cen(jc, 1);
              double dsqr = dx*dx + dy*dy;

              den += Densityfrac(jc, m) / dsqr;
              ++nnm;
            }
          }
          AvgMatDensity_row(m) = nnm > 0 ? den / nnm : 0.;
        }
        else
        {
           AvgMatDensity_row(m) = 0.;
        }
      }
    }

    timer.record();
    data_checker.check(MatDensity_average.getMap()->data());
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::neighbor_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_submap_templated, act_perf);

  SLIC_INFO("Average Material Density            compute time is "
            << act_perf << " secs\n");

}

//    Average material density over neighborhood of each cell
//      MultiMat - Dense - Iterator
void average_density_over_nbr_cell_dom_full_mm_iter(MultiMat& mm, Robey_data& data)
{
  SLIC_INFO("-- Average Density over Neighbors, Cell-Dominant Full Matrix,"
    << " Multimat Iterator --");
  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::CELL_DOM);
  SLIC_ASSERT(mm.getFieldSparsityLayout(0) == SparsityLayout::DENSE);

  const std::vector<double> & cen = data.cen;
  const std::vector<int> & nnbrs = data.nnbrs;
  const std::vector<int> & nbrs = data.nbrs;

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double> & Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double> & Volfrac = mm.get2dField<double>("Volfrac");

  std::vector<double> MatDensity_average(ncells * nmats, 0);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    for (auto& v : MatDensity_average)
      v = 0.0;

    timer.start();

    for (int ic = 0; ic < ncells; ++ic)
    {
      double xc[2];
      xc[0] = cen[ic * 2]; xc[1] = cen[ic * 2 + 1];
      int nn = nnbrs[ic];
      int cnbrs[8];
      double dsqr[8];

      for (int n = 0; n < nn; ++n)
        cnbrs[n] = nbrs[ic * 8 + n];

      for (int n = 0; n < nn; ++n)
      {
        dsqr[n] = 0.0;
        for (int d = 0; d < 2; ++d) //original condition was  d < 1 ???
        {
          double ddist = (xc[d] - cen[cnbrs[n] * 2 + d]);
          dsqr[n] += ddist * ddist;
        }
      }

      auto VolfracIter = Volfrac.begin(ic);
      auto VolfracIterEnd = Volfrac.end(ic);
      for( ;VolfracIter != VolfracIterEnd; ++VolfracIter)
      {
        int m = VolfracIter.index();
        if (*VolfracIter > 0.0)  //this check is not needed in sparse layout.
        {
          int nnm = 0;         // number of nbrs with this material
          for (int n = 0; n < nn; ++n)
          {
            int jc = cnbrs[n];
            auto VolfracIterJcm = Volfrac.begin(jc)+m;
            if (*(VolfracIterJcm) > 0.0)
            {
              auto DensityIterJcm = Densityfrac.begin(jc) + m;
              MatDensity_average[ic * nmats + m] += *DensityIterJcm / dsqr[n];
              ++nnm;
            }
          }
          if (nnm > 0)
            MatDensity_average[ic * nmats + m] /= nnm;
          else
            SLIC_ASSERT(MatDensity_average[ic * nmats + m] == 0.0);
        }
        else
        {
          MatDensity_average[ic * nmats + m] = 0.0;
        }
      }
    }

    timer.record();
    data_checker.check(MatDensity_average);
  } 

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::neighbor_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_iter, act_perf);
  SLIC_INFO("Average Material Density            compute time is "
    << act_perf << " secs\n");

}



//    Average material density over neighborhood of each cell
//      MultiMat - Compact - Submap
void average_density_over_nbr_cell_dom_compact_mm_submap(MultiMat& mm, Robey_data& data)
{
  SLIC_INFO("-- Average Density over Neighbors, Cell-Dominant Compact,"
    << " Multimat Submap --");
  mm.convertLayoutToCellDominant();
  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::CELL_DOM);
  SLIC_ASSERT(mm.getFieldSparsityLayout(0) == SparsityLayout::SPARSE);

  const auto& neighbors = data.slam_neighbors;
  const auto& cen = data.slam_centroids;

  int ncells = mm.getNumberOfCells();
  MultiMat::Field2D<double> & Densityfrac = mm.get2dField<double>("Densityfrac");
  //MultiMat::Field2D<double> & Volfrac = mm.get2dField<double>("Volfrac");
  auto MatDensity_average = mm.get2dField<double>("MatDensityAverage");

  //std::vector<double> MatDensity_average(ncells * nmats, 0);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    for (auto& v : MatDensity_average)
      v = 0.0;

    timer.start();

    for (int ic = 0; ic < ncells; ++ic)
    {
      auto nbrs = neighbors[ic];
      double xc[2] = { cen(ic,0), cen(ic,1) };

      auto Densityfrac_row = Densityfrac(ic);
      auto AvgMatDensity_row = MatDensity_average(ic);
      for (int k = 0; k < Densityfrac_row.size(); ++k)
      {
        const int m = Densityfrac_row.index(k);

        double den = 0.;
        int nnm = 0;         // number of nbrs with this material
        for (int n : nbrs.positions())
        {
          int jc = nbrs[n];
          const auto* val = Densityfrac.findValue(jc, m);
          if (val != nullptr)
          {
            double dx = xc[0] - cen(jc, 0);
            double dy = xc[1] - cen(jc, 1);
            double dsqr = dx*dx + dy * dy;

            den += *val / dsqr;
            ++nnm;
          }
        }

        AvgMatDensity_row(k) = nnm > 0 ? den / nnm : 0.;
      }
    }

    timer.record();
    data_checker.check(MatDensity_average.getMap()->data());
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::neighbor_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_submap, act_perf);
  SLIC_INFO("Average Material Density            compute time is "
    << act_perf << " secs\n");

}

template<typename BSet>
void average_density_over_nbr_cell_dom_compact_mm_submap(MultiMat& mm, Robey_data& data)
{
  SLIC_INFO("-- Average Density over Neighbors, Cell-Dominant Compact, "
      << " Multimat Submap -- templated on BSet type");
  mm.convertLayoutToCellDominant();
  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
      << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::CELL_DOM);
  SLIC_ASSERT(mm.getFieldSparsityLayout(0) == SparsityLayout::SPARSE);

  int ncells = mm.getNumberOfCells();

  //auto Densityfrac = mm.getTemplated2DField<double, DataLayout::CELL_DOM, BSet>("Densityfrac");
  //auto MatDensity_average = mm.getTemplated2DField<double, DataLayout::CELL_DOM, BSet>("MatDensityAverage");

  auto Densityfrac = mm.get2dFieldAsSlamBivarMap<double, BSet>("Densityfrac");
  auto MatDensity_average = mm.get2dFieldAsSlamBivarMap<double, BSet>("MatDensityAverage");

  // Get the slam relations and maps for the mesh data
  const auto& neighbors = data.slam_neighbors;
  const auto& cen = data.slam_centroids;

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    MatDensity_average.clear();

    timer.start();

    for (int ic = 0; ic < ncells; ++ic)
    {
      // Get the set of neighbors
      auto nbrs = neighbors[ic];

      // compute the squared distances to neighbors
      double xc[2] = { cen(ic,0), cen(ic,1) };

      auto Densityfrac_row = Densityfrac(ic);
      auto AvgMatDensity_row = MatDensity_average(ic);
      for (int k = 0; k < Densityfrac_row.size(); ++k)
      {
        const int m = Densityfrac_row.index(k);

        double den = 0.;
        int nnm = 0;         // number of nbrs with this material
        for (int n : nbrs.positions())
        {
          int jc = nbrs[n];
          const auto* val = Densityfrac.findValue(jc, m);
          if ( val != nullptr)
          {
            double dx = xc[0] - cen(jc, 0);
            double dy = xc[1] - cen(jc, 1);
            double dsqr = dx*dx + dy * dy;

            den += *val / dsqr;
            ++nnm;
          }
        }
                
        AvgMatDensity_row(k) = nnm > 0 ? den / nnm : 0.;
      }
    }

    timer.record();
    data_checker.check(MatDensity_average.getMap()->data());
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::neighbor_density, mm.getFieldDataLayout(0),
       mm.getFieldSparsityLayout(0), Result_Store::mm_submap_templated, act_perf);
  SLIC_INFO("Average Material Density            compute time is "
     << act_perf << " secs\n");
}

template<typename BSet, typename FieldTypeT>
void average_density_over_nbr_cell_dom_compact_mm_submap(MultiMat& mm, Robey_data& data)
{
  SLIC_INFO("-- Average Density over Neighbors, Cell-Dominant Compact, "
    << " Multimat Submap -- templated on BSet type");
  mm.convertLayoutToCellDominant();
  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::CELL_DOM);
  SLIC_ASSERT(mm.getFieldSparsityLayout(0) == SparsityLayout::SPARSE);

  int ncells = mm.getNumberOfCells();

  //pick which get field function to use
  GetFieldFunType<FieldTypeT> getFieldFun = getfieldfun<FieldTypeT, BSet>();

  FieldTypeT Densityfrac = (mm.*getFieldFun)("Densityfrac");
  FieldTypeT MatDensity_average = (mm.*getFieldFun)("MatDensityAverage");

  //auto Densityfrac = mm.getTemplated2DField<double, DataLayout::CELL_DOM, BSet>("Densityfrac");
  //auto MatDensity_average = mm.getTemplated2DField<double, DataLayout::CELL_DOM, BSet>("MatDensityAverage");

  //auto Densityfrac = mm.get2dFieldAsSlamBivarMap<double, BSet>("Densityfrac");
  //auto MatDensity_average = mm.get2dFieldAsSlamBivarMap<double, BSet>("MatDensityAverage");

  // Get the slam relations and maps for the mesh data
  const auto& neighbors = data.slam_neighbors;
  const auto& cen = data.slam_centroids;

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    MatDensity_average.clear();

    timer.start();

    for (int ic = 0; ic < ncells; ++ic)
    {
      // Get the set of neighbors
      auto nbrs = neighbors[ic];

      // compute the squared distances to neighbors
      double xc[2] = { cen(ic,0), cen(ic,1) };

      auto Densityfrac_row = Densityfrac(ic);
      auto AvgMatDensity_row = MatDensity_average(ic);
      for (int k = 0; k < Densityfrac_row.size(); ++k)
      {
        const int m = Densityfrac_row.index(k);

        double den = 0.;
        int nnm = 0;         // number of nbrs with this material
        for (int n : nbrs.positions())
        {
          int jc = nbrs[n];
          const auto* val = Densityfrac.findValue(jc, m);
          if (val != nullptr)
          {
            double dx = xc[0] - cen(jc, 0);
            double dy = xc[1] - cen(jc, 1);
            double dsqr = dx * dx + dy * dy;

            den += *val / dsqr;
            ++nnm;
          }
        }

        AvgMatDensity_row(k) = nnm > 0 ? den / nnm : 0.;
      }
    }

    timer.record();
    data_checker.check(MatDensity_average.getMap()->data());
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::neighbor_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_submap_templated, act_perf);
  SLIC_INFO("Average Material Density            compute time is "
    << act_perf << " secs\n");
}


//    Average material density over neighborhood of each cell
//      MultiMat - Compact - IndexArray
void average_density_over_nbr_cell_dom_compact_mm_idxarray(MultiMat& mm, Robey_data& data)
{
  SLIC_INFO("-- Average Density over Neighbors, Cell-Dominant Compact,"
    << " Multimat Index Array --");
  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::CELL_DOM);
  SLIC_ASSERT(mm.getFieldSparsityLayout(0) == SparsityLayout::SPARSE);

  const std::vector<double> & cen = data.cen;
  const std::vector<int> & nnbrs = data.nnbrs;
  const std::vector<int> & nbrs = data.nbrs;

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double> & Densityfrac = mm.get2dField<double>("Densityfrac");
  //MultiMat::Field2D<double> & Volfrac = mm.get2dField<double>("Volfrac");
  
  std::vector<double> MatDensity_average(ncells * nmats, 0);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    for (auto& v : MatDensity_average)
      v = 0.0;

    timer.start();

    for (int ic = 0; ic < ncells; ++ic)
    {
      double xc[2];
      xc[0] = cen[ic * 2]; xc[1] = cen[ic * 2 + 1];
      int nn = nnbrs[ic];
      int cnbrs[8];
      double dsqr[8];

      for (int n = 0; n < nn; ++n)
        cnbrs[n] = nbrs[ic * 8 + n];

      for (int n = 0; n < nn; ++n)
      {
        dsqr[n] = 0.0;
        for (int d = 0; d < 2; ++d) //original condition was  d < 1 ???
        {
          double ddist = (xc[d] - cen[cnbrs[n] * 2 + d]);
          dsqr[n] += ddist * ddist;
        }
      }

      auto IndexSet = mm.getIndexingSetOfCell(ic, mm.getFieldSparsityLayout(0));
      auto MatId = mm.getMatInCell(ic);
      
      for (int k = 0; k < IndexSet.size(); ++k)
      {
        int m = MatId[k];
        int nnm = 0;         // number of nbrs with this material
        for (int n = 0; n < nn; ++n)
        {
          int jc = cnbrs[n];
          auto MatIdJc = mm.getMatInCell(jc);
          for (int jj = 0; jj < MatIdJc.size(); ++jj)
          {
            if (MatIdJc[jj] == m)
            {
              auto IndexSetJc = mm.getIndexingSetOfCell(jc, mm.getFieldSparsityLayout(0));
              MatDensity_average[ic * nmats + m] += Densityfrac[IndexSetJc[jj]] / dsqr[n];
              ++nnm;
              break;
            }
          }
        }
        if (nnm > 0)
          MatDensity_average[ic * nmats + m] /= nnm;
        else
          SLIC_ASSERT(MatDensity_average[ic * nmats + m] == 0.0);

      }
    }

    timer.record();
    data_checker.check(MatDensity_average);
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::neighbor_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_idxarray, act_perf);
  SLIC_INFO("Average Material Density            compute time is "
    << act_perf << " secs\n");

}


//    Average material density over neighborhood of each cell
//      MultiMat - Compact - Iter
void average_density_over_nbr_cell_dom_compact_mm_iter(MultiMat& mm, Robey_data& data)
{
  SLIC_INFO("-- Average Density over Neighbors, Cell-Dominant Compact,"
    << " Multimat Iter --");
  mm.convertLayoutToCellDominant();
  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::CELL_DOM);
  SLIC_ASSERT(mm.getFieldSparsityLayout(0) == SparsityLayout::SPARSE);

  const std::vector<double> & cen = data.cen;
  const std::vector<int> & nnbrs = data.nnbrs;
  const std::vector<int> & nbrs = data.nbrs;

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double> & Densityfrac = mm.get2dField<double>("Densityfrac");
  
  std::vector<double> MatDensity_average(ncells * nmats, 0);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    for (auto& v : MatDensity_average)
      v = 0.0;

    timer.start();

    for (int ic = 0; ic < ncells; ++ic)
    {
      double xc[2];
      xc[0] = cen[ic * 2]; xc[1] = cen[ic * 2 + 1];
      int nn = nnbrs[ic];
      int cnbrs[8];
      double dsqr[8];

      for (int n = 0; n < nn; ++n)
        cnbrs[n] = nbrs[ic * 8 + n];

      for (int n = 0; n < nn; ++n)
      {
        dsqr[n] = 0.0;
        for (int d = 0; d < 2; ++d) //original condition was  d < 1 ???
        {
          double ddist = (xc[d] - cen[cnbrs[n] * 2 + d]);
          dsqr[n] += ddist * ddist;
        }
      }

      auto DensityfracIter = Densityfrac.begin(ic);
      auto DensityfracIterEnd = Densityfrac.end(ic);
      for (; DensityfracIter != DensityfracIterEnd ; ++DensityfracIter)
      {
        int m = DensityfracIter.index();

        int nnm = 0;         // number of nbrs with this material
        for (int n = 0; n < nn; ++n)
        {
          int jc = cnbrs[n];
          auto DensityfracIterJc = Densityfrac.begin(jc);
          auto DensityfracIterJcEnd = Densityfrac.end(jc);
          for (; DensityfracIterJc != DensityfracIterJcEnd; ++DensityfracIterJc)
          {
            if (DensityfracIterJc.index() == m)
            {
              MatDensity_average[ic * nmats + m] += *DensityfracIterJc / dsqr[n];
              ++nnm;
              break;
            }
          }
        }
        if (nnm > 0)
          MatDensity_average[ic * nmats + m] /= nnm;
        else
          SLIC_ASSERT(MatDensity_average[ic * nmats + m] == 0.0);

      }
    }

    timer.record();
    data_checker.check(MatDensity_average);
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::neighbor_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_iter, act_perf);
  SLIC_INFO("Average Material Density            compute time is "
    << act_perf << " secs\n");

}



////////////////////// Avg Neighbor - Material Dominant /////////////////////////

//    Average density over neighbor cells
//      Robey's - Material-Dominant Full Matrix
void average_density_over_nbr_mat_dom_full(Robey_data& data)
{
  int ncells = data.ncells;
  int nmats = data.nmats;
  std::vector<double>& Volfrac = data.Volfrac;
  std::vector<double>& Densityfrac = data.Densityfrac;
  const std::vector<double>& cen = data.cen;
  const std::vector<int>& nnbrs = data.nnbrs;
  const std::vector<int>& nbrs = data.nbrs;

  SLIC_INFO("-- Average Density over Neighbors, Material-Dominant Full Matrix Array Access --");

  std::vector<double> MatDensity_average(ncells * nmats, 0.0);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    for (auto& v : MatDensity_average)
      v = 0.0;

    timer.start();

    for (int m = 0; m < nmats; ++m) 
    {
      for (int ic = 0; ic < ncells; ++ic) 
      {
        if (Volfrac[m * ncells + ic] > 0.0) 
        {
          double xc[2] = { cen[ic * 2], cen[ic * 2 + 1] };
          int nn = nnbrs[ic];
          const int* cnbrs = &(nbrs[ic * 8]);

          int nnm = 0;         // number of nbrs with this material
          double den = 0.0;
          for (int n = 0; n < nn; ++n) {
            int jc = cnbrs[n];
            if (Volfrac[m*ncells+jc] > 0.0) 
            {
              double dx = xc[0] - cen[jc*2 + 0];
              double dy = xc[1] - cen[jc*2 + 1];
              double dsqr = dx*dx + dy*dy;

              den += Densityfrac[m * ncells + jc] / dsqr;
              ++nnm;
            }
          }
          if (nnm > 0)
            MatDensity_average[m * ncells + ic] = den / nnm;
          else
            SLIC_ASSERT(MatDensity_average[m * ncells + ic] == 0.0);
        }
        else {
          MatDensity_average[m * ncells + ic] = 0.0;
        }
      }
    }

    timer.record();
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::neighbor_density, DataLayout::MAT_DOM,
    SparsityLayout::DENSE, Result_Store::method_csr, act_perf);
  SLIC_INFO("Average Material Density            compute time is "
    << act_perf << " secs\n");

}


//    Average density over neighbor cells
//      Robey's - Material-Dominant Compact Matrix
void average_density_over_nbr_mat_dom_compact(Robey_data& data)
{
  int ncells = data.ncells;
  int nmats = data.nmats;
  //std::vector<double>& Volfrac = data.Volfrac_sparse;
  std::vector<double>& Densityfrac = data.Densityfrac_sparse;
  const std::vector<double>& cen = data.cen;
  const std::vector<int>& nnbrs = data.nnbrs;
  const std::vector<int>& nbrs = data.nbrs;
  std::vector<int>& begin_idx = data.begin_idx;
  std::vector<int>& cell_id = data.col_idx;

  SLIC_INFO("-- Average Density over Neighbors, Material-Dominant Compact Array Access --");

  std::vector<double> MatDensity_average(begin_idx.back(), 0.0);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    for (auto& v : MatDensity_average)
      v = 0.0;

    timer.start();

    for (int m = 0; m < nmats; ++m) {
      for (int ii = begin_idx[m]; ii < begin_idx[m + 1]; ++ii)
      {
        int ci = cell_id[ii];

        double xc[2] = { cen[ci * 2], cen[ci * 2 + 1] };
        int nn = nnbrs[ci];
        const int* cnbrs = &(nbrs[ci * 8]);
        
        int nnm = 0;         // number of nbrs with this material
        int begin_idx_m = begin_idx[m];
        double den = 0.0;
        for (int n = 0; n < nn; ++n) {
          int C_j = cnbrs[n];
          int c_j = data.dense2sparse_idx[m * ncells + C_j];
          if (c_j >= 0) { //the neighbor cell does have this material

            double dx = xc[0] - cen[C_j * 2 + 0];
            double dy = xc[1] - cen[C_j * 2 + 1];
            double dsqr = dx*dx + dy*dy;

            den += Densityfrac[begin_idx_m + c_j] / dsqr;
            ++nnm;
          }
        }
        if (nnm > 0)
          MatDensity_average[ii] = den / nnm;
        else
          SLIC_ASSERT(MatDensity_average[ii] == 0.0);
      }
    }

    timer.record();
    data_checker.check(MatDensity_average);
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::neighbor_density, DataLayout::MAT_DOM,
    SparsityLayout::SPARSE, Result_Store::method_csr, act_perf);
  SLIC_INFO("Average Material Density            compute time is "
    << act_perf << " secs\n");

}


//    Average density over neighbor cells
//      MultiMat - full - Direct Access
void average_density_over_nbr_mat_dom_full_mm_direct(MultiMat& mm, Robey_data& data)
{
  SLIC_INFO("-- Average Density over Neighbors, Material-Dominant, using MultiMat Direct Access --");

  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::MAT_DOM);
  SLIC_ASSERT(mm.getFieldSparsityLayout(0) == SparsityLayout::DENSE);

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double>& Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double>& Volfrac = mm.get2dField<double>("Volfrac");

  const std::vector<double>& cen = data.cen;
  const std::vector<int>& nnbrs = data.nnbrs;
  const std::vector<int>& nbrs = data.nbrs;

  std::vector<double> MatDensity_average(ncells * nmats, 0.0);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    for (auto& v : MatDensity_average)
      v = 0.0;

    timer.start();

    for (int m = 0; m < nmats; ++m) {
      for (int ic = 0; ic < ncells; ++ic) {
        if (*Volfrac.findValue(m, ic) > 0.0) {
          double xc[2];
          xc[0] = cen[ic * 2]; xc[1] = cen[ic * 2 + 1];
          int nn = nnbrs[ic];
          int cnbrs[8];
          double dsqr[8];
          for (int n = 0; n < nn; ++n) cnbrs[n] = nbrs[ic * 8 + n];
          for (int n = 0; n < nn; ++n) {
            dsqr[n] = 0.0;
            for (int d = 0; d < 2; ++d) { //????
              double ddist = (xc[d] - cen[cnbrs[n] * 2 + d]);
              dsqr[n] += ddist * ddist;
            }
          }

          int nnm = 0;         // number of nbrs with this material
          for (int n = 0; n < nn; ++n) {
            int jc = cnbrs[n];
            if (*Volfrac.findValue(m, jc) > 0.0) {
              MatDensity_average[m * ncells + ic] += *Densityfrac.findValue(m,jc) / dsqr[n];
              ++nnm;
            }
          }
          if (nnm > 0)
            MatDensity_average[m * ncells + ic] /= nnm;
          else
            SLIC_ASSERT(MatDensity_average[m * ncells + ic] == 0.0);
        }
        else {
          MatDensity_average[m * ncells + ic] = 0.0;
        }
      }
    }

    timer.record();
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::neighbor_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_direct, act_perf);
  SLIC_INFO("Average Material Density            compute time is "
    << act_perf << " secs\n");

}

//    Average density over neighbor cells
//      MultiMat - full - Submap
void average_density_over_nbr_mat_dom_full_mm_submap(MultiMat& mm, Robey_data& data)
{
  SLIC_INFO("-- Average Density over Neighbors, Material-Dominant, using MultiMat Submap --");

  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::MAT_DOM);
  SLIC_ASSERT(mm.getFieldSparsityLayout(0) == SparsityLayout::DENSE);

  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double> & Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double> & Volfrac = mm.get2dField<double>("Volfrac");
  auto MatDensity_average = mm.get2dField<double>("MatDensityAverage");

  // Get the slam relations and maps for the mesh data
  const auto& neighbors = data.slam_neighbors;
  const auto& cen = data.slam_centroids;

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    for (auto& v : MatDensity_average)
      v = 0.0;

    timer.start();

    for (int m = 0; m < nmats; ++m) 
    {
      auto Densityfrac_row = Densityfrac(m);
      auto Volfrac_row = Volfrac(m);
      auto AvgMatDensity_row = MatDensity_average(m);

      for (int ic = 0; ic < Volfrac_row.size(); ++ic) 
      {
        if (Volfrac_row(ic) > 0.0) 
        {
          auto nbrs = neighbors[ic];
          double xc[2] = { cen(ic,0), cen(ic,1) };


          int nnm = 0;         // number of nbrs with this material
          double den = 0.0;
          for (int n : nbrs.positions()) 
          {
            int jc = nbrs[n];
            if (Volfrac_row(jc) > 0.0) {

              double dx = xc[0] - cen(jc, 0);
              double dy = xc[1] - cen(jc, 1);
              double dsqr = dx*dx + dy*dy;

              den += Densityfrac_row(jc) / dsqr;
              ++nnm;
            }
          }
          if (nnm > 0)
            AvgMatDensity_row(ic) = den / nnm;
          else
            SLIC_ASSERT(AvgMatDensity_row(ic) == 0.0);
        }
        else {
          AvgMatDensity_row(ic) = 0.0;
        }
      }
    }

    timer.record();
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::neighbor_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_submap, act_perf);
  SLIC_INFO("Average Material Density            compute time is "
    << act_perf << " secs\n");

}

template<typename BSet>
void average_density_over_nbr_mat_dom_full_mm_submap(MultiMat& mm, Robey_data& data)
{
  SLIC_INFO("-- Average Density over Neighbors, Material-Dominant, using MultiMat Submap -- templated on BSet");

  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
      << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::MAT_DOM);
  SLIC_ASSERT(mm.getFieldSparsityLayout(0) == SparsityLayout::DENSE);

  int nmats = mm.getNumberOfMaterials();

  //auto Densityfrac = mm.getTemplated2DField<double, DataLayout::MAT_DOM, BSet>("Densityfrac");
  //auto Volfrac = mm.getTemplated2DField<double, DataLayout::MAT_DOM, BSet>("Volfrac");
  //auto MatDensity_average = mm.getTemplated2DField<double, DataLayout::MAT_DOM, BSet>("MatDensityAverage");

  auto Densityfrac = mm.get2dFieldAsSlamBivarMap<double,BSet>("Densityfrac");
  auto Volfrac = mm.get2dFieldAsSlamBivarMap<double,BSet>("Volfrac");
  auto MatDensity_average = mm.get2dFieldAsSlamBivarMap<double, BSet>("MatDensityAverage");

  // Get the slam relations and maps for the mesh data
  const auto& neighbors = data.slam_neighbors;
  const auto& cen = data.slam_centroids;


  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    MatDensity_average.clear();

    timer.start();

    for (int m = 0; m < nmats; ++m)
    {
      auto Densityfrac_row = Densityfrac(m);
      auto Volfrac_row = Volfrac(m);
      auto AvgMatDensity_row = MatDensity_average(m);

      for (int ic = 0; ic < Volfrac_row.size(); ++ic)
      {
        if (Volfrac_row(ic) > 0.0) 
        {
          double xc[2] = { cen(ic,0), cen(ic,1) };

          double den = 0.;
          int nnm = 0;         // number of nbrs with this material

          auto nbrs = neighbors[ic]; 
          for (int n : nbrs.positions())
          {
            const int jc = nbrs[n];
            if (Volfrac_row(jc) > 0.0) 
            {
              const double dx = xc[0] - cen(jc, 0);
              const double dy = xc[1] - cen(jc, 1);
              const double dsqr = dx*dx + dy * dy;

              den += Densityfrac_row(jc) / dsqr;
              ++nnm;
            }
          }

          AvgMatDensity_row(ic) = nnm > 0 ? den / nnm : 0.;
        }
        else 
        {
          AvgMatDensity_row(ic) = 0.;
        }
      }
    }

    timer.record();
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::neighbor_density, mm.getFieldDataLayout(0),
      mm.getFieldSparsityLayout(0), Result_Store::mm_submap_templated, act_perf);
  SLIC_INFO("Average Material Density            compute time is "
      << act_perf << " secs\n");
}


template<typename BSet, typename FieldTypeT>
void average_density_over_nbr_mat_dom_full_mm_submap(MultiMat& mm, Robey_data& data)
{
  SLIC_INFO("-- Average Density over Neighbors, Material-Dominant, using MultiMat Submap -- templated on BSet");

  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::MAT_DOM);
  SLIC_ASSERT(mm.getFieldSparsityLayout(0) == SparsityLayout::DENSE);

  int nmats = mm.getNumberOfMaterials();

  //pick which get field function to use
  GetFieldFunType<FieldTypeT> getFieldFun = getfieldfun<FieldTypeT, BSet>();

  FieldTypeT Densityfrac = (mm.*getFieldFun)("Densityfrac");
  FieldTypeT Volfrac = (mm.*getFieldFun)("Volfrac");
  FieldTypeT MatDensity_average = (mm.*getFieldFun)("MatDensityAverage");

  // Get the slam relations and maps for the mesh data
  const auto& neighbors = data.slam_neighbors;
  const auto& cen = data.slam_centroids;


  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    MatDensity_average.clear();

    timer.start();

    for (int m = 0; m < nmats; ++m)
    {
      auto Densityfrac_row = Densityfrac(m);
      auto Volfrac_row = Volfrac(m);
      auto AvgMatDensity_row = MatDensity_average(m);

      for (int ic = 0; ic < Volfrac_row.size(); ++ic)
      {
        if (Volfrac_row(ic) > 0.0)
        {
          double xc[2] = { cen(ic,0), cen(ic,1) };

          double den = 0.;
          int nnm = 0;         // number of nbrs with this material

          auto nbrs = neighbors[ic];
          for (int n : nbrs.positions())
          {
            const int jc = nbrs[n];
            if (Volfrac_row(jc) > 0.0)
            {
              const double dx = xc[0] - cen(jc, 0);
              const double dy = xc[1] - cen(jc, 1);
              const double dsqr = dx * dx + dy * dy;

              den += Densityfrac_row(jc) / dsqr;
              ++nnm;
            }
          }

          AvgMatDensity_row(ic) = nnm > 0 ? den / nnm : 0.;
        }
        else
        {
          AvgMatDensity_row(ic) = 0.;
        }
      }
    }

    timer.record();
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::neighbor_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_submap_templated, act_perf);
  SLIC_INFO("Average Material Density            compute time is "
    << act_perf << " secs\n");
}


//    Average density over neighbor cells
//      MultiMat - full - Iter
void average_density_over_nbr_mat_dom_full_mm_iter(MultiMat& mm, Robey_data& data)
{
  SLIC_INFO("-- Average Density over Neighbors, Material-Dominant, using MultiMat Iter --");

  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::MAT_DOM);
  SLIC_ASSERT(mm.getFieldSparsityLayout(0) == SparsityLayout::DENSE);

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double> & Densityfrac = mm.get2dField<double>("Densityfrac");
  MultiMat::Field2D<double> & Volfrac = mm.get2dField<double>("Volfrac");

  const std::vector<double> & cen = data.cen;
  const std::vector<int> & nnbrs = data.nnbrs;
  const std::vector<int> & nbrs = data.nbrs;

  std::vector<double> MatDensity_average(ncells * nmats, 0.0);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    for (auto& v : MatDensity_average)
      v = 0.0;

    timer.start();

    for (int m = 0; m < nmats; ++m) {
      auto VolfracIter = Volfrac.begin(m);
      auto VolfracIterEnd = Volfrac.end(m);
      auto DensityIterBegin = Densityfrac.begin(m);
      auto VolfracIterBegin = Volfrac.begin(m);
      for (; VolfracIter != VolfracIterEnd; ++VolfracIter) {
        int ic = VolfracIter.index();
        if (*VolfracIter > 0.0) {
          double xc[2];
          xc[0] = cen[ic * 2]; xc[1] = cen[ic * 2 + 1];
          int nn = nnbrs[ic];
          int cnbrs[8];
          double dsqr[8];
          for (int n = 0; n < nn; ++n) cnbrs[n] = nbrs[ic * 8 + n];
          for (int n = 0; n < nn; ++n) {
            dsqr[n] = 0.0;
            for (int d = 0; d < 2; ++d) { //????
              double ddist = (xc[d] - cen[cnbrs[n] * 2 + d]);
              dsqr[n] += ddist * ddist;
            }
          }

          int nnm = 0;         // number of nbrs with this material

          for (int n = 0; n < nn; ++n) {
            int jc = cnbrs[n];
            
            if (*(VolfracIterBegin + jc) > 0.0) {
              MatDensity_average[m * ncells + ic] += *(DensityIterBegin + jc) / dsqr[n];
              ++nnm;
            }
          }
          if (nnm > 0)
            MatDensity_average[m * ncells + ic] /= nnm;
          else
            SLIC_ASSERT(MatDensity_average[m * ncells + ic] == 0.0);
        }
        else {
          MatDensity_average[m * ncells + ic] = 0.0;
        }
      }
    }

    timer.record();
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::neighbor_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_iter, act_perf);
  SLIC_INFO("Average Material Density            compute time is "
    << act_perf << " secs\n");

}


//    Average density over neighbor cells
//      MultiMat - Compact - Submap
void average_density_over_nbr_mat_dom_compact_mm_submap(MultiMat& mm, Robey_data& data)
{
  SLIC_INFO("-- Average Density over Neighbors, Material-Dominant Compact, using MultiMat Submap --");

  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::MAT_DOM);
  SLIC_ASSERT(mm.getFieldSparsityLayout(0) == SparsityLayout::SPARSE);

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double>& Densityfrac = mm.get2dField<double>("Densityfrac");

  //std::vector<double> MatDensity_average(Densityfrac.totalSize(), 0.0);
  MultiMat::Field2D<double>& MatDensity_average = mm.get2dField<double>("MatDensityAverage");

  // Get the slam relations and maps for the mesh data
  const auto& neighbors = data.slam_neighbors;
  const auto& cen = data.slam_centroids;

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    MatDensity_average.clear();
    
    timer.start();

    for (int m = 0; m < nmats; ++m) 
    {
      auto Densityfrac_row = Densityfrac(m);
      auto AvgMatDensity_row = MatDensity_average(m);
      for (int k = 0; k < Densityfrac_row.size(); ++k)
      {
        int ci = Densityfrac_row.index(k);

        auto nbrs = neighbors[ci];
        double xc[2] = { cen(ci,0), cen(ci,1) };

        int nnm = 0;         // number of nbrs with this material
        double den = 0.0;
        for (int n : nbrs.positions()) {
          int C_j = nbrs[n];
          int k_j = data.dense2sparse_idx[m * ncells + C_j];
          if (k_j >= 0) { //the neighbor cell does have this material
            double dx = xc[0] - cen[C_j * 2 + 0];
            double dy = xc[1] - cen[C_j * 2 + 1];
            double dsqr = dx*dx + dy*dy;

            den += Densityfrac_row(k_j) / dsqr;
            ++nnm;
          }
        }
        if (nnm > 0)
          AvgMatDensity_row(k) = den / nnm;
        else
          SLIC_ASSERT(AvgMatDensity_row(k) == 0.0);
      }
    }

    timer.record();
    data_checker.check(MatDensity_average.getMap()->data());
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::neighbor_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_submap, act_perf);
  SLIC_INFO("Average Material Density            compute time is "
    << act_perf << " secs\n");

}

template<typename BSet>
void average_density_over_nbr_mat_dom_compact_mm_submap(MultiMat& mm, Robey_data& data)
{
  SLIC_INFO("-- Average Density over Neighbors, Material-Dominant Compact, using MultiMat Submap -- templated on BSet");

  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
      << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::MAT_DOM);
  SLIC_ASSERT(mm.getFieldSparsityLayout(0) == SparsityLayout::SPARSE);

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  
  //auto Densityfrac = mm.getTemplated2DField<double, DataLayout::MAT_DOM, BSet>("Densityfrac");
  //auto MatDensity_average = mm.getTemplated2DField<double, DataLayout::MAT_DOM, BSet>("MatDensityAverage");

  auto Densityfrac = mm.get2dFieldAsSlamBivarMap<double,BSet>("Densityfrac");
  auto MatDensity_average = mm.get2dFieldAsSlamBivarMap<double, BSet>("MatDensityAverage");

  // Get the slam relations and maps for the mesh data
  const auto& neighbors = data.slam_neighbors;
  const auto& cen = data.slam_centroids;

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    MatDensity_average.clear();

    timer.start();

    for (int m = 0; m < nmats; ++m)
    {
      auto Densityfrac_row = Densityfrac(m);
      auto AvgMatDensity_row = MatDensity_average(m);

      for (int k = 0; k < Densityfrac_row.size(); ++k)
      {
        int ci = Densityfrac_row.index(k);

        double xc[2] = { cen(ci,0), cen(ci,1) };

        // Get the set of neighbors
        auto nbrs = neighbors[ci];

        double den = 0.;
        int nnm = 0;         // number of nbrs with this material
        for (int n : nbrs.positions() ) 
        {
          int C_j = nbrs[n];
          int k_j = data.dense2sparse_idx[m * ncells + C_j];
          if (k_j >= 0)  //the neighbor cell has this material
          {
            const double dx = xc[0] - cen(C_j, 0);
            const double dy = xc[1] - cen(C_j, 1);
            const double dsqr = dx*dx + dy * dy;

            den += Densityfrac_row(k_j) / dsqr;
            ++nnm;
          }
        }

        AvgMatDensity_row(k) = nnm > 0 ? den / nnm : 0;
      }
    }

    timer.record();
    data_checker.check(MatDensity_average.getMap()->data());
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::neighbor_density, mm.getFieldDataLayout(0),
     mm.getFieldSparsityLayout(0), Result_Store::mm_submap_templated, act_perf);
  SLIC_INFO("Average Material Density            compute time is "
      << act_perf << " secs\n");
}


template<typename BSet, typename FieldTypeT>
void average_density_over_nbr_mat_dom_compact_mm_submap(MultiMat& mm, Robey_data& data)
{
  SLIC_INFO("-- Average Density over Neighbors, Material-Dominant Compact, using MultiMat Submap -- templated on BSet");

  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::MAT_DOM);
  SLIC_ASSERT(mm.getFieldSparsityLayout(0) == SparsityLayout::SPARSE);

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();

  //pick which get field function to use
  GetFieldFunType<FieldTypeT> getFieldFun = getfieldfun<FieldTypeT, BSet>();

  FieldTypeT Densityfrac = (mm.*getFieldFun)("Densityfrac");
  FieldTypeT MatDensity_average = (mm.*getFieldFun)("MatDensityAverage");

  // Get the slam relations and maps for the mesh data
  const auto& neighbors = data.slam_neighbors;
  const auto& cen = data.slam_centroids;

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    MatDensity_average.clear();

    timer.start();

    for (int m = 0; m < nmats; ++m)
    {
      auto Densityfrac_row = Densityfrac(m);
      auto AvgMatDensity_row = MatDensity_average(m);

      for (int k = 0; k < Densityfrac_row.size(); ++k)
      {
        int ci = Densityfrac_row.index(k);

        double xc[2] = { cen(ci,0), cen(ci,1) };

        // Get the set of neighbors
        auto nbrs = neighbors[ci];

        double den = 0.;
        int nnm = 0;         // number of nbrs with this material
        for (int n : nbrs.positions())
        {
          int C_j = nbrs[n];
          int k_j = data.dense2sparse_idx[m * ncells + C_j];
          if (k_j >= 0)  //the neighbor cell has this material
          {
            const double dx = xc[0] - cen(C_j, 0);
            const double dy = xc[1] - cen(C_j, 1);
            const double dsqr = dx * dx + dy * dy;

            den += Densityfrac_row(k_j) / dsqr;
            ++nnm;
          }
        }

        AvgMatDensity_row(k) = nnm > 0 ? den / nnm : 0;
      }
    }

    timer.record();
    data_checker.check(MatDensity_average.getMap()->data());
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::neighbor_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_submap_templated, act_perf);
  SLIC_INFO("Average Material Density            compute time is "
    << act_perf << " secs\n");
}

//    Average density over neighbor cells
//      MultiMat - Compact - Index Array
void average_density_over_nbr_mat_dom_compact_mm_indexarray(MultiMat& mm, Robey_data& data)
{
  SLIC_INFO("-- Average Density over Neighbors, Material-Dominant Compact, using MultiMat Index Array--");

  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::MAT_DOM);
  SLIC_ASSERT(mm.getFieldSparsityLayout(0) == SparsityLayout::SPARSE);

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double> & Densityfrac = mm.get2dField<double>("Densityfrac");
  //MultiMat::Field2D<double> & Volfrac = mm.get2dField<double>("Volfrac");

  const std::vector<double> & cen = data.cen;
  const std::vector<int> & nnbrs = data.nnbrs;
  const std::vector<int> & nbrs = data.nbrs;

  std::vector<double> MatDensity_average(ncells * nmats, 0.0);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    for (auto& v : MatDensity_average)
      v = 0.0;

    timer.start();

    for (int m = 0; m < nmats; ++m)
    {
      auto IndexSet = mm.getIndexingSetOfMat(m, SparsityLayout::SPARSE);
      auto CellIdArr = mm.getCellContainingMat(m);

      for (int k = 0; k < IndexSet.size(); ++k)
      {
        int ci = CellIdArr[k];

        double xc[2];
        xc[0] = cen[ci * 2]; xc[1] = cen[ci * 2 + 1];
        int nn = nnbrs[ci];
        int cnbrs[9];
        double dsqr[8];
        for (int n = 0; n < nn; ++n) cnbrs[n] = nbrs[ci * 8 + n];
        for (int n = 0; n < nn; ++n) {
          dsqr[n] = 0.0;
          for (int d = 0; d < 2; ++d) { //???
            double ddist = (xc[d] - cen[cnbrs[n] * 2 + d]);
            dsqr[n] += ddist * ddist;
          }
        }

        int nnm = 0;         // number of nbrs with this material
        for (int n = 0; n < nn; ++n) {
          int C_j = cnbrs[n];
          int k_j = data.dense2sparse_idx[m * ncells + C_j];
          if (k_j >= 0) { //the neighbor cell does have this material
            MatDensity_average[m * ncells + ci] += Densityfrac[IndexSet[k_j]] / dsqr[n];
            ++nnm;
          }
        }
        if (nnm > 0)
          MatDensity_average[m * ncells + ci] /= nnm;
        else
          SLIC_ASSERT(MatDensity_average[m * ncells + ci] == 0.0);
      }
    }

    timer.record();
    data_checker.check(MatDensity_average);
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::neighbor_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_idxarray, act_perf);
  SLIC_INFO("Average Material Density            compute time is "
    << act_perf << " secs\n");

}

//    Average density over neighbor cells
//      MultiMat - Compact - Iter
void average_density_over_nbr_mat_dom_compact_mm_iter(MultiMat& mm, Robey_data& data)
{
  SLIC_INFO("-- Average Density over Neighbors, Material-Dominant compact, using MultiMat Iter --");

  SLIC_INFO("MultiMat layout: " << mm.getFieldDataLayoutAsString(0) << " & "
    << mm.getFieldSparsityLayoutAsString(0));
  SLIC_ASSERT(mm.getFieldDataLayout(0) == DataLayout::MAT_DOM);
  SLIC_ASSERT(mm.getFieldSparsityLayout(0) == SparsityLayout::SPARSE);

  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();
  MultiMat::Field2D<double>& Densityfrac = mm.get2dField<double>("Densityfrac");

  const std::vector<double>& cen = data.cen;
  const std::vector<int>& nnbrs = data.nnbrs;
  const std::vector<int>& nbrs = data.nbrs;

  std::vector<double> MatDensity_average(ncells * nmats, 0.0);

  timer.reset();

  for (int iter = 0; iter < ITERMAX; ++iter)
  {
    for (auto& v : MatDensity_average)
      v = 0.0;

    timer.start();

    for (int m = 0; m < nmats; ++m)
    {
      auto DensityfracIter = Densityfrac.begin(m);
      auto DensityfracIterBegin = Densityfrac.begin(m);
      auto DensityfracIterEnd = Densityfrac.end(m);
      for (;DensityfracIter!= DensityfracIterEnd; ++DensityfracIter)
      {
        int ci = DensityfracIter.index();

        double xc[2];
        xc[0] = cen[ci * 2]; xc[1] = cen[ci * 2 + 1];
        int nn = nnbrs[ci];
        int cnbrs[9];
        double dsqr[8];
        for (int n = 0; n < nn; ++n) cnbrs[n] = nbrs[ci * 8 + n];
        for (int n = 0; n < nn; ++n) {
          dsqr[n] = 0.0;
          for (int d = 0; d < 2; ++d) { //???
            double ddist = (xc[d] - cen[cnbrs[n] * 2 + d]);
            dsqr[n] += ddist * ddist;
          }
        }

        int nnm = 0;         // number of nbrs with this material
        for (int n = 0; n < nn; ++n) {
          int C_j = cnbrs[n];
          int k_j = data.dense2sparse_idx[m * ncells + C_j];
          if (k_j >= 0) { //the neighbor cell does have this material
            MatDensity_average[m * ncells + ci] += *(DensityfracIterBegin + k_j) / dsqr[n];
            ++nnm;
          }
        }
        if (nnm > 0)
          MatDensity_average[m * ncells + ci] /= nnm;
        else
          SLIC_ASSERT(MatDensity_average[m * ncells + ci] == 0.0);
      }
    }

    timer.record();
    data_checker.check(MatDensity_average);
  }

  double act_perf = timer.get_median();
  result_store.add_result(Result_Store::neighbor_density, mm.getFieldDataLayout(0),
    mm.getFieldSparsityLayout(0), Result_Store::mm_iter, act_perf);
  SLIC_INFO("Average Material Density            compute time is "
    << act_perf << " secs\n");

}






////////////////////// Main /////////////////////////

int main(int argc, char** argv)
{
  axom::slic::UnitTestLogger logger;

  if (argc != 1 && argc != 2 && argc != 3)
  {
    SLIC_WARNING("Usage: ./multimat_calculate_ex [<volfrac_data>]");
    SLIC_WARNING("  or   ./multimat_calculate_ex <ncells> <nmats>");
    return 1;
  }

  std::string fileName = "";
  int ncells_to_gen = 1600;
  int nmats_to_gen = 10;

  if (argc == 2)
  {
    fileName = std::string(argv[1]);
  }
  else if (argc == 3)
  {
    ncells_to_gen = atoi(argv[1]);
    nmats_to_gen = atoi(argv[2]);
  }

  Robey_data data(fileName, ncells_to_gen, nmats_to_gen);
  data.set_up_cell_dom_data();

  result_store.init(&data);

  //Set-up the multimat class
  MultiMat mm(DataLayout::CELL_DOM, SparsityLayout::DENSE);
  mm.setNumberOfMaterials(data.nmats);
  mm.setNumberOfCells(data.ncells);
  mm.setCellMatRel(data.Volfrac_bool, DataLayout::CELL_DOM);

  //Setting field data in terms of slam
  mm.addField<>("Densityfrac", FieldMapping::PER_CELL_MAT, DataLayout::CELL_DOM, SparsityLayout::DENSE, &data.Densityfrac[0]);
  mm.addField<>("Vol", FieldMapping::PER_CELL, DataLayout::CELL_DOM, SparsityLayout::DENSE, &data.Vol[0]);
  mm.addField<>("Volfrac", FieldMapping::PER_CELL_MAT, DataLayout::CELL_DOM, SparsityLayout::DENSE, &data.Volfrac[0]);
  mm.addField<>("Tempfrac", FieldMapping::PER_CELL_MAT, DataLayout::CELL_DOM, SparsityLayout::DENSE, &data.Temperaturefrac[0]);
  mm.addField<>("Pressurefrac", FieldMapping::PER_CELL_MAT, DataLayout::CELL_DOM, SparsityLayout::DENSE, &data.Pressurefrac[0]);
  mm.addField<>("nmatconsts", FieldMapping::PER_MAT, DataLayout::CELL_DOM, SparsityLayout::DENSE, &data.nmatconsts[0]);
  mm.addField<>("MatDensityAverage", FieldMapping::PER_CELL_MAT, DataLayout::CELL_DOM, SparsityLayout::DENSE, &data.Pressurefrac[0]);


  //printself and check
  mm.isValid(true);

  using SetType = slam::RangeSet<>;
  //using ProductSetType = slam::ProductSet<SetType, SetType>;
  using ProductSetType = MultiMat::ProductSetType; 

  using RelType = MultiMat::SparseRelationType;
  //using RelationSetType = slam::RelationSet<RelType, SetType, SetType>;
  using RelationSetType = MultiMat::RelationSetType;

  // test out field by field change
  mm.convertFieldToCellDom(1);
  std::cout << "Layout after CD-convert: " << mm.getFieldDataLayoutAsString(1) << std::endl;

  auto df = mm.get2dField<double>("Densityfrac");
  std::cout << "first set size: " << df.firstSetSize() << std::endl;

  mm.convertFieldToMatDom(1);
  std::cout << "Layout after MD-convert: " << mm.getFieldDataLayoutAsString(1) << std::endl;

  df = mm.get2dField<double>("Densityfrac");
  std::cout << "first set size: " << df.firstSetSize() << std::endl;

  mm.convertFieldToCellDom(1);
  std::cout << "Layout after CD-convert: " << mm.getFieldDataLayoutAsString(1) << std::endl;

  df = mm.get2dField<double>("Densityfrac");
  std::cout << "first set size: " << df.firstSetSize() << std::endl;

  mm.convertFieldToMatDom(1);
  std::cout << "Layout after MD-convert: " << mm.getFieldDataLayoutAsString(1) << std::endl;

  df = mm.get2dField<double>("Densityfrac");
  std::cout << "first set size: " << df.firstSetSize() << std::endl;


  SLIC_INFO("**********************************************************");
  SLIC_INFO("* ");
  SLIC_INFO("* Average Density");
  SLIC_INFO("* ");

  //Run the Full layout, cell dom
  SLIC_INFO("*************** Full Layout - Cell Dominant **************");
  data_checker.reset();
  mm.convertLayoutToCellDominant();
  mm.convertLayoutToDense();

  average_density_cell_dom_full(data);
  average_density_cell_dom_mm_template(mm);
  average_density_cell_dom_mm_direct(mm);

  average_density_cell_dom_mm_submap(mm);
#ifdef run_slam_bivarmap
  average_density_cell_dom_mm_submap<ProductSetType>(mm);
  average_density_cell_dom_mm_submap<ProductSetType, BiVarMapT<ProductSetType>>(mm);
#endif
  average_density_cell_dom_mm_submap<ProductSetType, Field2DT<ProductSetType>>(mm);
  average_density_cell_dom_mm_submap<ProductSetType, Field2DTempT<DataLayout::CELL_DOM, ProductSetType>>(mm);
  average_density_cell_dom_mm_idxarray(mm);

  average_density_cell_dom_mm_iter(mm);
  //average_density_cell_dom_mm_flatiter(mm);

  //Run the Compact layout, cell dom
  SLIC_INFO("*************** Compact Layout - Cell Dominant **************");
  average_density_cell_dom_compact(data);
  mm.convertLayoutToSparse();
  average_density_cell_dom_mm_idxarray(mm);

  average_density_cell_dom_mm_submap(mm);
#ifdef run_slam_bivarmap
  average_density_cell_dom_mm_submap<RelationSetType>(mm);
  average_density_cell_dom_mm_submap<RelationSetType, BiVarMapT<RelationSetType>>(mm);
#endif
  average_density_cell_dom_mm_submap<RelationSetType, Field2DT<RelationSetType>>(mm);
  average_density_cell_dom_mm_submap<RelationSetType, Field2DTempT<DataLayout::CELL_DOM, RelationSetType>>(mm);

  average_density_cell_dom_mm_iter(mm);
  //average_density_cell_dom_mm_flatiter(mm);

  //return 0;

  //Run the Full layout, material dom
  SLIC_INFO("*************** Full Layout - Material Dominant **************");
  data_checker.reset();

  data.set_up_mat_dom_data();
  mm.convertLayoutToDense();
  mm.convertLayoutToMaterialDominant();

  average_density_mat_dom_full(data);

  average_density_mat_dom_mm_direct(mm);
  average_density_mat_dom_mm_submap(mm);
#ifdef run_slam_bivarmap
  average_density_mat_dom_mm_submap<ProductSetType>(mm);
  average_density_mat_dom_mm_submap<ProductSetType, BiVarMapT<ProductSetType>>(mm);
#endif
  average_density_mat_dom_mm_submap<ProductSetType, Field2DT<ProductSetType>>(mm);
  average_density_mat_dom_mm_submap<ProductSetType, Field2DTempT<DataLayout::MAT_DOM, ProductSetType>>(mm);
  average_density_mat_dom_mm_iter(mm);
  //average_density_mat_dom_mm_flatiter(mm);
  
  SLIC_INFO("*************** Compact Layout - Material Dominant **************");
  mm.convertLayoutToSparse();

  average_density_mat_dom_compact(data); 
  average_density_mat_dom_mm_idxarray(mm);
  average_density_mat_dom_mm_submap(mm);
#ifdef run_slam_bivarmap
  average_density_mat_dom_mm_submap<RelationSetType>(mm);
  average_density_mat_dom_mm_submap<RelationSetType, BiVarMapT<RelationSetType>>(mm);
#endif
  average_density_mat_dom_mm_submap<RelationSetType, Field2DT<RelationSetType>>(mm);
  average_density_mat_dom_mm_submap<RelationSetType, Field2DTempT<DataLayout::MAT_DOM, RelationSetType>>(mm);

  average_density_mat_dom_mm_iter(mm);
  //average_density_mat_dom_mm_flatiter(mm);

  SLIC_INFO("**********************************************************");
  SLIC_INFO("* ");
  SLIC_INFO("* Average Density over Neighbor");
  SLIC_INFO("* ");

  SLIC_INFO("**************** Full Layout - Cell-Dominant ******************");
  data.set_up_cell_dom_data();
  data_checker.reset();
  mm.convertLayoutToCellDominant();
  mm.convertLayoutToDense();
  average_density_over_nbr_cell_dom_full(data);
  average_density_over_nbr_cell_dom_full_mm_direct(mm, data);
  average_density_over_nbr_cell_dom_full_mm_submap(mm, data);
#ifdef run_slam_bivarmap
  average_density_over_nbr_cell_dom_full_mm_submap<ProductSetType>(mm, data);
  average_density_over_nbr_cell_dom_full_mm_submap<ProductSetType, BiVarMapT<ProductSetType>>(mm, data);
#endif
  average_density_over_nbr_cell_dom_full_mm_submap<ProductSetType, Field2DT<ProductSetType>>(mm, data);
  average_density_over_nbr_cell_dom_full_mm_submap<ProductSetType, Field2DTempT<DataLayout::CELL_DOM, ProductSetType>>(mm, data);
  average_density_over_nbr_cell_dom_full_mm_iter(mm, data);

  SLIC_INFO("**************** Compact Layout - Cell-Dominant ******************");
  mm.convertLayoutToSparse();
  data_checker.reset();
  average_density_over_nbr_cell_dom_compact(data);
  //average_density_over_nbr_cell_dom_compact_mm_idxarray(mm, data);
  average_density_over_nbr_cell_dom_compact_mm_submap(mm, data);
#ifdef run_slam_bivarmap
  average_density_over_nbr_cell_dom_compact_mm_submap<RelationSetType>(mm, data);
  average_density_over_nbr_cell_dom_compact_mm_submap<RelationSetType, BiVarMapT<RelationSetType>>(mm, data);
#endif
  average_density_over_nbr_cell_dom_compact_mm_submap<RelationSetType, Field2DT<RelationSetType>>(mm, data);
  average_density_over_nbr_cell_dom_compact_mm_submap<RelationSetType, Field2DTempT<DataLayout::CELL_DOM, RelationSetType>>(mm, data);
  //average_density_over_nbr_cell_dom_compact_mm_iter(mm, data);

  SLIC_INFO("**************** Full Layout - Material-Dominant ******************");
  data.set_up_mat_dom_data();
  data_checker.reset();
  mm.convertLayoutToMaterialDominant();
  mm.convertLayoutToDense();

  average_density_over_nbr_mat_dom_full(data);
  average_density_over_nbr_mat_dom_full_mm_direct(mm, data);
  average_density_over_nbr_mat_dom_full_mm_submap(mm, data);
#ifdef run_slam_bivarmap
  average_density_over_nbr_mat_dom_full_mm_submap<ProductSetType>(mm, data);
  average_density_over_nbr_mat_dom_full_mm_submap<ProductSetType, BiVarMapT<ProductSetType>>(mm, data);
#endif
  average_density_over_nbr_mat_dom_full_mm_submap<ProductSetType, Field2DT<ProductSetType>>(mm, data);
  average_density_over_nbr_mat_dom_full_mm_submap<ProductSetType, Field2DTempT<DataLayout::MAT_DOM, ProductSetType>>(mm, data);
  //average_density_over_nbr_mat_dom_full_mm_iter(mm, data);

  SLIC_INFO("**************** Compact Layout - Material-Dominant ******************");
  mm.convertLayoutToSparse();
  data_checker.reset();
  average_density_over_nbr_mat_dom_compact(data);

  //average_density_over_nbr_mat_dom_compact_mm_indexarray(mm, data);
  average_density_over_nbr_mat_dom_compact_mm_submap(mm, data);
#ifdef run_slam_bivarmap
  average_density_over_nbr_mat_dom_compact_mm_submap<RelationSetType>(mm, data);
  average_density_over_nbr_mat_dom_compact_mm_submap<RelationSetType, BiVarMapT<RelationSetType>>(mm, data);
#endif
  average_density_over_nbr_mat_dom_compact_mm_submap<RelationSetType, Field2DT<RelationSetType>>(mm, data);
  average_density_over_nbr_mat_dom_compact_mm_submap<RelationSetType, Field2DTempT<DataLayout::MAT_DOM, RelationSetType>>(mm, data);
  //average_density_over_nbr_mat_dom_compact_mm_iter(mm, data);

  SLIC_INFO("**********************************************************");
  SLIC_INFO("* ");
  SLIC_INFO("* Calculate Pressure Using Ideal Gas Law");
  SLIC_INFO("* ");
  
  SLIC_INFO("**************** Full Layout - Cell-Dominant ******************");

  data.set_up_cell_dom_data();
  data_checker.reset();
  mm.convertLayoutToCellDominant();
  mm.convertLayoutToDense();
  calculate_pressure_cell_dom_full(data);
  calculate_pressure_cell_dom_full_mm_direct(mm);
  calculate_pressure_cell_dom_mm_submap(mm);
#ifdef run_slam_bivarmap
  calculate_pressure_cell_dom_mm_submap<ProductSetType>(mm);
  calculate_pressure_cell_dom_mm_submap<ProductSetType, BiVarMapT<ProductSetType>>(mm);
#endif
  calculate_pressure_cell_dom_mm_submap<ProductSetType, Field2DT<ProductSetType>>(mm);
  calculate_pressure_cell_dom_mm_submap<ProductSetType, Field2DTempT<DataLayout::CELL_DOM, ProductSetType>>(mm);
  calculate_pressure_cell_dom_full_mm_iter(mm);
  //calculate_pressure_cell_dom_full_mm_flatiter(mm);

  SLIC_INFO("**************** Compact Layout - Cell-Dominant ******************");
  mm.convertLayoutToSparse();
  data_checker.reset();
  calculate_pressure_cell_dom_compact(data);

  calculate_pressure_cell_dom_mm_submap(mm);
#ifdef run_slam_bivarmap
  calculate_pressure_cell_dom_mm_submap<RelationSetType>(mm);
  calculate_pressure_cell_dom_mm_submap<RelationSetType, BiVarMapT<RelationSetType>>(mm);
#endif
  calculate_pressure_cell_dom_mm_submap<RelationSetType, Field2DT<RelationSetType>>(mm);
  calculate_pressure_cell_dom_mm_submap<RelationSetType, Field2DTempT<DataLayout::CELL_DOM, RelationSetType>>(mm);

  SLIC_INFO("**************** Full Layout - Material-Dominant ******************");
  data.set_up_mat_dom_data();
  data_checker.reset();
  mm.convertLayoutToDense();
  mm.convertLayoutToMaterialDominant();

  calculate_pressure_mat_dom_full(data);
  calculate_pressure_mat_dom_full_mm_direct(mm);
  calculate_pressure_mat_dom_mm_submap(mm);
#ifdef run_slam_bivarmap
  calculate_pressure_mat_dom_mm_submap<ProductSetType>(mm);
  calculate_pressure_mat_dom_mm_submap<ProductSetType, BiVarMapT<ProductSetType>>(mm);
#endif
  calculate_pressure_mat_dom_mm_submap<ProductSetType, Field2DT<ProductSetType>>(mm);
  calculate_pressure_mat_dom_mm_submap<ProductSetType, Field2DTempT<DataLayout::MAT_DOM, ProductSetType>>(mm);
  calculate_pressure_mat_dom_full_mm_iter(mm);
  //calculate_pressure_mat_dom_full_mm_flatiter(mm);

  SLIC_INFO("**************** Compact Layout - Material-Dominant ******************");
  mm.convertLayoutToSparse();
  calculate_pressure_mat_dom_compact(data);

  calculate_pressure_mat_dom_mm_submap(mm);
#ifdef run_slam_bivarmap
  calculate_pressure_mat_dom_mm_submap<RelationSetType>(mm);
  calculate_pressure_mat_dom_mm_submap<RelationSetType, BiVarMapT<RelationSetType>>(mm);
#endif
  calculate_pressure_mat_dom_mm_submap<RelationSetType, Field2DT<RelationSetType>>(mm);
  calculate_pressure_mat_dom_mm_submap<RelationSetType, Field2DTempT<DataLayout::MAT_DOM, RelationSetType>>(mm);

  /*
  average_density_cell_dom_with_if(data);
  average_density_cell_dom_with_if_mm(mm);
  */


  SLIC_INFO("*");
  SLIC_INFO("*  Finished.");
  SLIC_INFO("**********************************************************");
  SLIC_INFO("");

  result_store.save_to_csv_file("test_save_result.csv");

  return 0;
}
