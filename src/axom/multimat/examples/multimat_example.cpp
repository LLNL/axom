
#include "multimat_example.hpp"

#include "helper.hpp"
#include <ctime>



using namespace std;
using namespace axom::multimat;

//using MM_doubleArrType = MultiMatTypedArray<double>;

#define ITERMAX 10  //define how many iterations to run test code

int method = 0;     // VF initialization: 0 - random, 1 - read volfrac.dat
float filled_fraction = -1;

void test_code();

std::clock_t start_time;

void start_timer() {
  start_time = std::clock();
}
double end_timer() {
  double duration = (std::clock() - start_time) / (double)CLOCKS_PER_SEC;
  return duration;
}
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////    Average density with fractional densities
////    Copied from Robey's code for Cell-Dominant Full Matrix Data Structure line 257
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void average_density_cell_dom(int ncells, int nmats, std::vector<double>& Volfrac, std::vector<double>& Densityfrac, vector<double>& Vol) {
//  cout << "-- Averaging Density cell-dominant array-access --" << endl;
//  vector<double> Density_average(ncells);
//
//  double time_sum = 0;
//  for (int iter = 0; iter < ITERMAX; iter++) {
//    std::clock_t start;
//    double duration;
//    start = std::clock();
//
//    for (int ic = 0; ic < ncells; ic++) {
//      double density_ave = 0.0;
//      for (int m = 0; m < nmats; m++) {
//        density_ave += Densityfrac[ic*nmats + m] * Volfrac[ic*nmats + m];
//      }
//      Density_average[ic] = density_ave / Vol[ic];
//    }
//
//    duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
//    time_sum += duration;
//
//    //Do a check for calculation correctness
//    double sum = 0;
//    for (auto da : Density_average)
//      sum += da;
//    cout << "iter" << iter << "   sum of values: " << sum << endl;
//  }
//
//  double act_perf = time_sum / (double)ITERMAX;
//  printf("Average Density of mixed material cells    compute time is %lf secs\n", act_perf);
//
//}
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////    Average density with fractional densities
////    Same as the function above, but modified to use MultiMat class
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void average_density_cell_dom_mm(MultiMat& mm) {
//  cout << "-- Averaging Density cell-dominant mm-version --" << endl;
//  int ncells = mm.m_ncells;
//  int nmats = mm.m_nmats;
//  MM_doubleArrType* Densityfrac = MM_CAST_TO(double, mm.getFieldArray("Densityfrac"));
//  MM_doubleArrType* Volfrac = MM_CAST_TO(double, mm.getFieldArray("Volfrac"));
//  MM_doubleArrType* Vol = MM_CAST_TO(double, mm.getFieldArray("Vol"));
//
//  vector<double> Density_average(ncells);
//
//  double time_sum = 0;
//
//  for (int iter = 0; iter < ITERMAX; iter++) {
//    std::clock_t start;
//    double duration;
//    start = std::clock();
//
//    for (int ic = 0; ic < ncells; ic++) {
//      double density_ave = 0.0;
//      for (int m = 0; m < nmats; m++) {
//        density_ave += Densityfrac->getValue(ic, m) *  Volfrac->getValue(ic, m);
//      }
//      Density_average[ic] = density_ave / Vol->getValue(ic);
//    }
//
//    duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
//    time_sum += duration;
//
//    //Do a check for calculation correctness
//    double sum = 0;
//    for (auto da : Density_average)
//      sum += da;
//    cout << "iter" << iter << "   sum of values: " << sum << endl;
//  }
//  double act_perf = time_sum / (double)ITERMAX;
//  printf("Average Density of mixed material cells    compute time is %lf secs\n", act_perf);
//}
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////    Average density with fractional densities
////        Material-first loop
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void average_density_mat_dom_mm(MultiMat& mm) {
//  cout << "-- Averaging Density mat-dominant mm-version --" << endl;
//  int ncells = mm.m_ncells;
//  int nmats = mm.m_nmats;
//  MM_doubleArrType* Densityfrac = MM_CAST_TO(double, mm.getFieldArray("Densityfrac"));;
//  MM_doubleArrType* Volfrac = MM_CAST_TO(double, mm.getFieldArray("Volfrac"));
//  MM_doubleArrType* Vol = MM_CAST_TO(double, mm.getFieldArray("Vol"));
//
//  double time_sum = 0;
//
//  for (int iter = 0; iter < ITERMAX; iter++) {
//    vector<double> Density_average(ncells);
//
//    std::clock_t start;
//    double duration;
//    start = std::clock();
//
//    for (int m = 0; m < nmats; m++) {
//      for (int ic = 0; ic < ncells; ic++) {
//        Density_average[ic] += Densityfrac->getValue(ic, m) * Volfrac->getValue(ic, m);
//      }
//    }
//    for (int ic = 0; ic < ncells; ic++) {
//      Density_average[ic] /= Vol->getValue(ic);
//    }
//
//    duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
//    time_sum += duration;
//
//    //Do a check for calculation correctness
//    double sum = 0;
//    for (auto da : Density_average)
//      sum += da;
//    cout << "iter" << iter << "   sum of values: " << sum << endl;
//  }
//  double act_perf = time_sum / (double)ITERMAX;
//  printf("Average Density of mixed material cells    compute time is %lf secs\n", act_perf);
//}
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////    Average density with fractional densities
////    Same as the function above, but modified to use MultiMat class
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void average_density_cell_dom_iter_mm(MultiMat& mm) {
//  cout << "-- Averaging Density cell-dominant using iterator mm-version --" << endl;
//  //TODO when mm class has iterators
//}
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////    Average density with fractional densities with if test
////    Copied from Robey's code for Cell-Dominant Full Matrix Data Structure
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void average_density_cell_dom_with_if(int ncells, int nmats, std::vector<double>& Volfrac, std::vector<double>& Densityfrac, vector<double>& Vol)
//{
//  cout << "-- Averaging Density with if --" << endl;
//  vector<double> Density_average(ncells, 0);
//  double time_sum = 0;
//
//  for (int iter = 0; iter < ITERMAX; iter++) {
//    std::clock_t start;
//    double duration;
//    start = std::clock();
//
//    for (int ic = 0; ic < ncells; ic++) {
//      double density_ave = 0.0;
//      for (int m = 0; m < nmats; m++) {
//        if (Volfrac[ic*nmats + m] > 0.0) {
//          density_ave += Densityfrac[ic*nmats + m] * Volfrac[ic*nmats + m];
//        }
//      }
//      Density_average[ic] = density_ave / Vol[ic];
//    }
//
//    duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
//    time_sum += duration;
//
//    //Do a check for calculation correctness
//    double sum = 0;
//    for (auto da : Density_average)
//      sum += da;
//    cout << "iter" << iter << "   sum of values: " << sum << endl;
//  }
//
//  double act_perf = time_sum / (double)ITERMAX;
//  printf("Average Density of frac with if            compute time is %lf secs\n", act_perf);
//
//}
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////    Average density with fractional densities with if test
////          Same as the function above, but modified to use MultiMat class
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void average_density_cell_dom_with_if_mm(MultiMat& mm) {
//  cout << "-- Averaging Density with if, mm-version --" << endl;
//  int ncells = mm.m_ncells;
//  int nmats = mm.m_nmats;
//  MM_doubleArrType* Densityfrac = MM_CAST_TO(double, mm.getFieldArray("Densityfrac"));;
//  MM_doubleArrType* Volfrac = MM_CAST_TO(double, mm.getFieldArray("Volfrac"));
//  MM_doubleArrType* Vol = MM_CAST_TO(double, mm.getFieldArray("Vol"));
//
//  vector<double> Density_average(ncells, 0);
//
//  double time_sum = 0;
//
//  for (int iter = 0; iter < ITERMAX; iter++) {
//    std::clock_t start;
//    double duration;
//    start = std::clock();
//
//    for (int ic = 0; ic < ncells; ic++) {
//      double density_ave = 0.0;
//      for (int m = 0; m < nmats; m++) {
//        if (Volfrac->getValue(ic, m) > 0.0) {
//          density_ave += Densityfrac->getValue(ic, m) *  Volfrac->getValue(ic, m);
//        }
//      }
//      Density_average[ic] = density_ave / Vol->getValue(ic);
//    }
//
//    duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
//    time_sum += duration;
//
//    //Do a check for calculation correctness
//    double sum = 0;
//    for (auto da : Density_average)
//      sum += da;
//    cout << "iter" << iter << "   sum of values: " << sum << endl;
//  }
//
//  double act_perf = time_sum / (double)ITERMAX;
//  printf("Average Density of frac with if            compute time is %lf secs\n", act_perf);
//
//}
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////   Calculate pressure using ideal gas law
////    Copied from Robey's code for Cell-Dominant Full Matrix Data Structure
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void calculate_pressure(int ncells, int nmats, std::vector<double>& Volfrac, std::vector<double>& Densityfrac,
//  vector<double>& Temperaturefrac, vector<double>& nmatconsts)
//{
//  cout << "-- Calculating pressure --" << endl;
//  vector<double> Pressurefrac(ncells*nmats, 0);
//
//  double time_sum = 0;
//
//  for (int iter = 0; iter< ITERMAX; iter++) {
//    std::clock_t start;
//    double duration;
//    start = std::clock();
//
//    for (int ic = 0; ic < ncells; ic++) {
//      for (int m = 0; m < nmats; m++) {
//        if (Volfrac[ic*nmats + m] > 0.) {
//          Pressurefrac[ic*nmats + m] = (nmatconsts[m] * Densityfrac[ic*nmats + m] * Temperaturefrac[ic*nmats + m]) / (Volfrac[ic*nmats + m]);
//        }
//        else {
//          Pressurefrac[ic*nmats + m] = 0.0;
//        }
//      }
//    }
//
//    duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
//    time_sum += duration;
//
//    //Do a check for calculation correctness
//    double sum = 0;
//    for (auto da : Pressurefrac)
//      sum += da;
//    cout << "iter" << iter << "   sum of values: " << sum << endl;
//  }
//  double act_perf = time_sum / (double)ITERMAX;
//  printf("Pressure Calculation of mixed material cells with if compute time is %lf msecs\n", act_perf);
//
//}
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////   Calculate pressure using ideal gas law
////          Same as the function above, but modified to use MultiMat class
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void calculate_pressure_mm(MultiMat& mm)
//{
//  cout << "-- Calculating pressure, mm-version --" << endl;
//  int ncells = mm.m_ncells;
//  int nmats = mm.m_nmats;
//  MM_doubleArrType* Densityfrac = MM_CAST_TO(double, mm.getFieldArray("Densityfrac"));;
//  MM_doubleArrType* Volfrac = MM_CAST_TO(double, mm.getFieldArray("Volfrac"));
//  //  MM_doubleArrType* Vol = MM_CAST_TO(double, mm.getFieldArray("Vol"));
//  MM_doubleArrType* Temperaturefrac = MM_CAST_TO(double, mm.getFieldArray("Temperaturefrac"));
//
//  vector<double> Pressurefrac(ncells*nmats, 0);
//
//  double time_sum = 0;
//
//  for (int iter = 0; iter< ITERMAX; iter++) {
//    std::clock_t start;
//    double duration;
//    start = std::clock();
//
//    for (int ic = 0; ic < ncells; ic++) {
//      for (int m = 0; m < nmats; m++) {
//        if (Volfrac->getValue(ic, m) > 0.) {
//          Pressurefrac[ic*nmats + m] = ( /*nmatconsts[m]*/ 5.0 * Densityfrac->getValue(ic, m)
//            * Temperaturefrac->getValue(ic, m)) / Volfrac->getValue(ic, m);
//        }
//        else {
//          Pressurefrac[ic*nmats + m] = 0.0;
//        }
//      }
//    }
//
//    duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
//    time_sum += duration;
//
//    //Do a check for calculation correctness
//    double sum = 0;
//    for (auto da : Pressurefrac)
//      sum += da;
//    cout << "iter" << iter << "   sum of values: " << sum << endl;
//  }
//  double act_perf = time_sum / (double)ITERMAX;
//  printf("Pressure Calculation of mixed material cells with if compute time is %lf msecs\n", act_perf);
//
//}
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////    Average material density over neighborhood of each cell
////    Copied from Robey's code for Cell-Dominant Full Matrix Data Structure
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void average_material_density_over_cell_nbr(int ncells, int nmats, std::vector<double>& Volfrac,
//  std::vector<double>& Densityfrac, //vector<double>& Temperaturefrac, vector<double>& nmatconsts
//  const std::vector<double>& cen, const std::vector<int>& nnbrs, const std::vector<int>& nbrs)
//{
//  cout << "-- Calculating avg material density over cell neighbor --" << endl;
//
//  double time_sum = 0;
//
//  for (int iter = 0; iter< ITERMAX; iter++) {
//    vector<double> MatDensity_average(ncells*nmats, 0);
//
//    std::clock_t start;
//    double duration;
//    start = std::clock();
//
//    for (int ic = 0; ic < ncells; ic++) {
//      double xc[2];
//      xc[0] = cen[ic * 2]; xc[1] = cen[ic * 2 + 1];
//      int nn = nnbrs[ic];
//      int cnbrs[8];
//      double dsqr[8];
//
//      for (int n = 0; n < nn; n++)
//        cnbrs[n] = nbrs[ic * 8 + n];
//
//      for (int n = 0; n < nn; n++) {
//        dsqr[n] = 0.0;
//        for (int d = 0; d < 2; d++) { //1???
//          double ddist = (xc[d] - cen[cnbrs[n] * 2 + d]);
//          dsqr[n] += ddist * ddist;
//        }
//      }
//
//      for (int m = 0; m < nmats; m++) {
//        if (Volfrac[ic*nmats + m] > 0.0) {
//          int nnm = 0;         // number of nbrs with this material
//          for (int n = 0; n < nn; n++) {
//            int jc = cnbrs[n];
//            if (Volfrac[ic*nmats + m] > 0.0) {
//              MatDensity_average[ic*nmats + m] += Densityfrac[ic*nmats + m] / dsqr[n];
//              nnm++;
//            }
//          }
//          MatDensity_average[ic*nmats + m] /= nnm;
//        }
//        else {
//          MatDensity_average[ic*nmats + m] = 0.0;
//        }
//      }
//    }
//
//    duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
//    time_sum += duration;
//
//    //Do a check for calculation correctness
//    double sum = 0;
//    for (auto da : MatDensity_average)
//      sum += da;
//    cout << "iter" << iter << "   sum of values: " << sum << endl;
//  }
//
//  double act_perf = time_sum / (double)ITERMAX;
//
//  printf("Average Material Density            compute time is %lf msecs\n", act_perf);
//
//}
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////    Average material density over neighborhood of each cell
////    Same as the function above, but modified to use MultiMat class
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void average_material_density_over_cell_nbr_mm(MultiMat& mm, const std::vector<double>& cen, const std::vector<int>& nnbrs, const std::vector<int>& nbrs)
//{
//  cout << "-- Calculating avg material density over cell neighbor, Multimat version--" << endl;
//
//  int ncells = mm.m_ncells;
//  int nmats = mm.m_nmats;
//  MM_doubleArrType* Densityfrac = MM_CAST_TO(double, mm.getFieldArray("Densityfrac"));;
//  MM_doubleArrType* Volfrac = MM_CAST_TO(double, mm.getFieldArray("Volfrac"));
//  //MM_doubleArrType* Vol = MM_CAST_TO(double, mm.getFieldArray("Vol"));
//  //MM_doubleArrType* Temperaturefrac = MM_CAST_TO(double, mm.getFieldArray("Temperaturefrac"));
//
//  double time_sum = 0;
//
//  for (int iter = 0; iter< ITERMAX; iter++) {
//    vector<double> MatDensity_average(ncells*nmats, 0);
//
//    std::clock_t start;
//    double duration;
//    start = std::clock();
//
//    for (int ic = 0; ic < ncells; ic++) {
//      double xc[2];
//      xc[0] = cen[ic * 2]; xc[1] = cen[ic * 2 + 1];
//      int nn = nnbrs[ic];
//      int cnbrs[8];
//      double dsqr[8];
//
//      for (int n = 0; n < nn; n++)
//        cnbrs[n] = nbrs[ic * 8 + n];
//
//      for (int n = 0; n < nn; n++) {
//        dsqr[n] = 0.0;
//        for (int d = 0; d < 2; d++) { //1???
//          double ddist = (xc[d] - cen[cnbrs[n] * 2 + d]);
//          dsqr[n] += ddist * ddist;
//        }
//      }
//
//      for (int m = 0; m < nmats; m++) {
//        if (Volfrac->getValue(ic, m)> 0.0) {
//          int nnm = 0;         // number of nbrs with this material
//          for (int n = 0; n < nn; n++) {
//            int jc = cnbrs[n];
//            if (Volfrac->getValue(ic, m) > 0.0) {
//              MatDensity_average[ic*nmats + m] += Densityfrac->getValue(ic, m) / dsqr[n];
//              nnm++;
//            }
//          }
//          MatDensity_average[ic*nmats + m] /= nnm;
//        }
//        else {
//          MatDensity_average[ic*nmats + m] = 0.0;
//        }
//      }
//    }
//
//    duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
//    time_sum += duration;
//
//    //Do a check for calculation correctness
//    double sum = 0;
//    for (auto da : MatDensity_average)
//      sum += da;
//    cout << "iter" << iter << "   sum of values: " << sum << endl;
//  }
//
//  double act_perf = time_sum / (double)ITERMAX;
//
//  printf("Average Material Density            compute time is %lf msecs\n", act_perf);
//
//}
//



int main(int argc, char** argv)
{
  int ncells;
  int nmats;
  float filled_percentage;

  vector<double> Volfrac;
  vector<double> Vol;
  vector<double> Densityfrac;
  vector<double> Temperaturefrac;
  vector<double> Pressurefrac;
  vector<double> nmatconsts;

  if (method) {
    //read from file... large and takes a long time.
    read_vol_frac_matrix_file(ncells, nmats, Volfrac, filled_percentage);
  }
  else {
    //create random data
    //get_vol_frac_matrix_rand(ncells, nmats, Volfrac, filled_percentage); 
    get_vol_frac_matrix_rand(ncells, nmats, Volfrac, filled_percentage, 100, 50); //small version
  }
  make_other_field_data(ncells, nmats, Volfrac, Vol, Densityfrac, Temperaturefrac, Pressurefrac, nmatconsts);
  filled_fraction = filled_percentage / 100.0;

  // Some variables on neighbors
  float L_f = method ? 0.5 : 1.0;  // ave frac of nbrs containing material
  int nnbrs_ave = 8;  // nearly so; 4000 boundary cells in 1 million cells
            //                  // in 3D, nnbrs_ave would be 26
  int nnbrs_max = 8;
  // Build up list of neighbors for each cell
  // Assuming a 2D structured mesh, each cell will have a maximum of 8 nbrs
  vector<int> nnbrs(ncells);
  vector<int> nbrs(ncells*nnbrs_max);
  get_neighbors(ncells, nnbrs_max, nnbrs, nbrs);

  // Compute centroids of cells
  vector<double> cen(ncells * 2);
  get_centroids(ncells, cen);

  //Making data for SLAM cell to mat relation
  vector<bool> Volfrac_bool(ncells*nmats, false);
  int cellmatcount = 0;
  for (int i = 0; i< Volfrac.size(); i++)
    if (Volfrac[i] > 0) {
      Volfrac_bool[i] = true;
      cellmatcount++;
    }

  //Get the non-zero fields for SALM set-up
  std::vector<double> Densityfrac_sparse(cellmatcount);
  std::vector<double> Volfrac_sparse(cellmatcount);
  std::vector<double> Temperaturefrac_sparse(cellmatcount);
  int ii = 0;
  for (int ci = 0; ci < ncells; ci++) {
    for (int mi = 0; mi < nmats; mi++) {
      double v = Densityfrac[ci*nmats + mi];
      if (v != 0) {
        Densityfrac_sparse[ii] = v;
        Volfrac_sparse[ii] = Volfrac[ci*nmats + mi];
        Temperaturefrac_sparse[ii] = Temperaturefrac[ci*nmats + mi];
        ii++;
      }
    }
  }

  //Set-up the multimat class
  MultiMat mm;
  mm.setNumberOfMat(nmats);
  mm.setNumberOfCell(ncells);
  mm.setCellMatRel(Volfrac_bool);


  //Setting field data in terms of slam
  auto MMArr_densityfrac     = mm.newFieldArray<>("Densityfrac"    , PER_CELL_MAT, &Densityfrac_sparse[0]);
  auto MMArr_vol             = mm.newFieldArray<>("Vol"            , PER_CELL    , &Vol[0]);
  auto MMArr_volfrac         = mm.newFieldArray<>("Volfrac"        , PER_CELL_MAT, &Volfrac_sparse[0]);
  auto MMArr_temperaturefrac = mm.newFieldArray<>("Temperaturefrac", PER_CELL_MAT, &Temperaturefrac_sparse[0]);
  auto MMArr_pressurefrac    = mm.newFieldArray<>("Pressurefrac"   , PER_CELL_MAT, &Pressurefrac[0]);
  auto MMArr_nmatconsts      = mm.newFieldArray<>("nmatconsts"     , PER_MAT     , &nmatconsts[0]);

  //printself and check
  mm.printSelf();
  printf("IsValid: %d\n\n", mm.isValid());

  //data value check, to make sure the values in multimat class is set-up correctly
  //for (int ci = 0; ci<mm.m_ncells; ci++) {
  //  for (int mi = 0; mi < mm.m_nmats; mi++) {
  //    double v = MMArr_densityfrac->getValue(ci, mi);
  //    assert(v == Densityfrac[ci*mm.m_nmats + mi]);

  //    v = MMArr_volfrac->getValue(ci, mi);
  //    assert(v == Volfrac[ci*mm.m_nmats + mi]);

  //    v = MMArr_temperaturefrac->getValue(ci, mi);
  //    assert(v == Temperaturefrac[ci*mm.m_nmats + mi]);

  //    v = MMArr_pressurefrac->getValue(ci, mi);
  //    assert(v == Pressurefrac[ci*mm.m_nmats + mi]);
  //  }
  //  assert(MMArr_vol->getValue(ci) == Vol[ci]);
  //}
  //for (int mi = 0; mi < mm.m_nmats; mi++) {
  //  double v = MMArr_nmatconsts->getValue(mi);
  //  assert(v == nmatconsts[mi]);
  //}


  ////Run the examples
  //average_density_cell_dom(ncells, nmats, Volfrac, Densityfrac, Vol);
  //average_density_cell_dom_mm(mm);
  //average_density_mat_dom_mm(mm);

  //average_density_cell_dom_with_if(ncells, nmats, Volfrac, Densityfrac, Vol);
  //average_density_cell_dom_with_if_mm(mm);

  //calculate_pressure(ncells, nmats, Volfrac, Densityfrac, Temperaturefrac, nmatconsts);
  //calculate_pressure_mm(mm);

  //average_material_density_over_cell_nbr(ncells, nmats, Volfrac, Densityfrac, cen, nnbrs, nbrs);
  //average_material_density_over_cell_nbr_mm(mm, cen, nnbrs, nbrs);

  test_code();

  cout << "\nPress Enter to terminate...";
  cin.ignore();
}


void test_code() {

  //Set-up
  MultiMat mm;
  int nmats = 100;
  int ncells = 20000;
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
  std::vector<double> cell_arr1(ncells);
  double c_sum = 0;
  for (int i = 0; i < cell_arr1.size(); i++) {
    cell_arr1[i] = (double)i * 2.0;
    c_sum += cell_arr1[i];
  }
    
  std::vector<double> cellmat_arr(nfilled);
  double x_sum = 0;
  for (int i = 0; i < cellmat_arr.size(); i++) {
    cellmat_arr[i] = (double)i * 1.1;
    x_sum += cellmat_arr[i];
  }


  auto MMarr_cell = mm.newFieldArray<>("Cell Array", PER_CELL, &cell_arr1[0]);
  //MMarr_cell->setValue(0, 123);
  //c_sum += 123;

  auto MMarr_cellmat = mm.newFieldArray<>("CellMat Array", PER_CELL_MAT, &cellmat_arr[0]);

  double sum = 0;

  //Different accessing methods ...


  // ------- Accessing using MMArray ----------
  printf("\nAccess through MMArray object (dense access)\n");
  sum = 0;
  start_timer();
  {
    MultiMatArray* cellarr = mm.getFieldArray("Cell Array");
    for (int i = 0; i < mm.m_ncells; i++) {
      sum += cellarr->getValue<double>(i);                  //<----
    }
  }
  cout << end_timer() << "\t";
  assert(c_sum == sum);

  sum = 0;
  start_timer();
  {
    MultiMatArray* cellmatArr = mm.getFieldArray("CellMat Array");
    for (int i = 0; i < mm.m_ncells; i++) {
      for (int m = 0; m < mm.m_nmats; m++) {
        sum += cellmatArr->getValue<double>(i, m);          //<---- warning: extra for-loop in the call!
      }
    }
  }
  cout << end_timer() << "\t";
  assert(x_sum == sum);

    
  // ------- Accessing using MMArray typed----------
  printf("\nAccess through MMTypedArray (dense access)\n");
  sum = 0;
  start_timer();
  {
    MultiMatTypedArray<double>* cellarr = mm.getFieldArray<double>("Cell Array");
    for (int i = 0; i < mm.m_ncells; i++) {
      sum += cellarr->getValue(i);                  //<----
    }
  }
  cout << end_timer() << "\t";
  assert(c_sum == sum);

  sum = 0;
  start_timer();
  {
    MultiMatTypedArray<double>* cellmatArr = mm.getFieldArray<double>("CellMat Array");
    for (int i = 0; i < mm.m_ncells; i++) {
      for (int m = 0; m < mm.m_nmats; m++) {
        sum += cellmatArr->getValue(i, m);          //<---- warning: extra for-loop in the call!
      }
    }
  }
  cout << end_timer() << "\t";
  assert(x_sum == sum);

  
  // --------- returning SLAM map -----------
  printf("\nAccess from SLAM map\n");
  sum = 0;
  start_timer();
  {
    //get it from MultiMai with array name
    //const MultiMat::MapType<double>& map = mm.getMap<double>("Cell Array");
   
    //MultiMatArray* cellarr = mm.getFieldArray("Cell Array");

    //Get a MapBaseType and cast it
    //const MultiMatArray::MapBaseType* map1 = cellarr->getMap();
    //const MultiMatArray::MapType<double>& map = *dynamic_cast<const MultiMatArray::MapType<double>*>(map1);

    //from MMArray
    //const MultiMatArray::MapType<double>& map = cellarr->getMap<double>();

    MultiMatTypedArray<double>* cellarr = mm.getFieldArray<double>("Cell Array");
    const MultiMatArray::MapType<double>& map = cellarr->getMap();
    
    for (int i = 0; i < mm.m_ncells; i++) {
      sum += map[i];                               //<----
    }
  }
  cout << end_timer() << "\t";
  assert(c_sum == sum);
  
  //don't have a way to access cell-mat relation through map type
  

  // ------------ return index set --------------
  printf("\nAccess by Map with indexing set\n-\t");
  sum = 0;
  start_timer();
  {
    MultiMatTypedArray<double>* cellmatarr = mm.getFieldArray<double>("CellMat Array");
    auto map = cellmatarr->getMap();
    for (int i = 0; i < mm.m_ncells; i++)
    {
      MultiMat::RelationSet setOfMaterialsInThisCell = mm.getMatInCell(i); //the materials (by id) in this cell
      MultiMat::RangeSetType indexSet = cellmatarr->getIndexingSetOfCell(i); //the indices into the maps
      assert(setOfMaterialsInThisCell.size() == indexSet.size());
      for (int j = 0; j < indexSet.size(); j++)
      {
        int mat_id = setOfMaterialsInThisCell[j];
        double val = map[indexSet.at(j)];
        sum += val;
      }
    }
  }
  cout << end_timer() << "\t";
  assert(x_sum == sum);


  // -------------- return subuset map for typed--------------

  printf("\nGet subset map from TypedArray\n-\t");
  sum = 0;
  start_timer();
  {
    MultiMatTypedArray<double>* cellmatarr = mm.getFieldArray<double>("CellMat Array");
    for (int i = 0; i < mm.m_ncells; i++)
    {
      MultiMat::RelationSet setOfMaterialsInThisCell = mm.getMatInCell(i);
      const MultiMatArray::MapType<double>& ValueMap = cellmatarr->getSubsetMap(i); //a map (basically std::vector) of the field for all the materials present in this cell
      assert(setOfMaterialsInThisCell.size() == ValueMap.size());
      for (int j = 0; j < ValueMap.size(); j++)
      {
        int mat_id = setOfMaterialsInThisCell[j];
        sum += ValueMap[j];
      }
    }
  }
  cout << end_timer() << "\t";
  assert(x_sum == sum);


  // -------------- return subuset map --------------

  printf("\nGet subset map\n-\t");
  sum = 0;
  start_timer();
  {
    MultiMatArray* cellmatarr = mm.getFieldArray("CellMat Array");
    for (int i = 0; i < mm.m_ncells; i++)
    {
      MultiMat::RelationSet setOfMaterialsInThisCell = mm.getMatInCell(i);
      const MultiMatArray::MapType<double>& ValueMap = cellmatarr->getSubsetMap<double>(i); //a map (basically std::vector) of the field for all the materials present in this cell
      assert(setOfMaterialsInThisCell.size() == ValueMap.size());
      for (int j = 0; j < ValueMap.size(); j++)
      {
        int mat_id = setOfMaterialsInThisCell[j];
        sum += ValueMap[j];
      }
    }
  }
  cout << end_timer() << "\t";
  assert(x_sum == sum);


  // ---------- using iterator -------------
  //template for all the iterator functions, which can be clunky
  printf("\nWith iterators\n");
  sum = 0;
  start_timer();
  {
    MultiMatArray* mm_cellarr = mm.getFieldArray("Cell Array");
    for (MultiMatArray::iterator<double> a = mm_cellarr->begin<double>(); a != mm_cellarr->end<double>(); a++) {
      sum += *a;              //<----
    }
  }
  cout << end_timer() << "\t";
  assert(c_sum == sum);

  sum = 0;
  start_timer();
  {
    MultiMatArray* mm_matcellarr = mm.getFieldArray("CellMat Array");
    for (MultiMatArray::iterator<double> a = mm_matcellarr->begin<double>(); a != mm_matcellarr->end<double>(); a++) {
      sum += *a;            //<----
    }
  }
  cout << end_timer() << "\t";
  assert(x_sum == sum);



  // ---------- using iterator with MultiMatTypedArray -------------
  //template for all the iterator functions, which can be clunky
  printf("\nWith MultiMatTypedArray iterators\n");
  sum = 0;
  start_timer();
  {
    MultiMatTypedArray<double>* cellarr = mm.getFieldArray<double>("Cell Array");
    for (MultiMatTypedArray<double>::iterator a = cellarr->begin(); a != cellarr->end(); a++) {
      sum += *a;              //<----
    }
  }
  cout << end_timer() << "\t";
  assert(c_sum == sum);

  sum = 0;
  start_timer();
  {
    MultiMatArray* mm_matcellarr = mm.getFieldArray("CellMat Array");
    for (MultiMatArray::iterator<double> a = mm_matcellarr->begin<double>(); a != mm_matcellarr->end<double>(); a++) {
      sum += *a;            //<----
    }
  }
  cout << end_timer() << "\t";
  assert(x_sum == sum);


  cout << endl;
}