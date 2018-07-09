
#include "multimat/multimat.hpp"

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

  std::vector<double> Densityfrac_sparse;
  std::vector<double> Volfrac_sparse;
  std::vector<double> Temperaturefrac_sparse;
  std::vector<double> Pressurefrac_sparse;

  Robey_data(bool read_from_file_bool = false)
  {

    if (read_from_file_bool) {
      //read from file... large and takes a long time.
      read_vol_frac_matrix_file(ncells, nmats, Volfrac, filled_percentage);
    }
    else {
      //create random data
      //get_vol_frac_matrix_rand(ncells, nmats, Volfrac, filled_percentage); 
      get_vol_frac_matrix_rand(ncells, nmats, Volfrac, filled_percentage, 100, 50); //small version
    }
    make_other_field_data(ncells, nmats, Volfrac, Vol, Densityfrac, Temperaturefrac, Pressurefrac, nmatconsts);
    filled_fraction = filled_percentage / 100.0f;

    // Some variables on neighbors
    float L_f = read_from_file_bool ? 0.5 : 1.0;  // ave frac of nbrs containing material
    int nnbrs_ave = 8;  // nearly so; 4000 boundary cells in 1 million cells
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

//
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////    Average density with fractional densities
////    Copied from Robey's code for Cell-Dominant Full Matrix Data Structure line 257
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void average_density_cell_dom(Robey_data data){
//  int ncells = data.ncells;
//  int nmats = data.nmats;
//  std::vector<double>& Volfrac = data.Volfrac;
//  std::vector<double>& Densityfrac = data.Densityfrac;
//  std::vector<double>& Vol = data.Vol;
//
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
//  int ncells = mm.getNumberOfCells();
//  int nmats = mm.getNumberOfMaterials();
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
//  int ncells = mm.getNumberOfCells();
//  int nmats = mm.getNumberOfMaterials();
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
//void average_density_cell_dom_with_if(Robey_data data)
//{
//  int ncells = data.ncells;
//  int nmats = data.nmats;
//  std::vector<double>& Volfrac = data.Volfrac;
//  std::vector<double>& Densityfrac = data.Densityfrac;
//  vector<double>& Vol = data.Vol;
//
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
//  int ncells = mm.getNumberOfCells();
//  int nmats = mm.getNumberOfMaterials();
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
//void calculate_pressure(Robey_data data)
//{
//  int ncells = data.ncells;
//  int nmats = data.nmats;
//  std::vector<double>& Volfrac = data.Volfrac;
//  std::vector<double>& Densityfrac = data.Densityfrac;
//
//  vector<double>& Temperaturefrac = data.Temperaturefrac;
//  vector<double>& nmatconsts = data.nmatconsts;
//
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
//  int ncells = mm.getNumberOfCells();
//  int nmats = mm.getNumberOfMaterials();
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
//void average_material_density_over_cell_nbr(Robey_data data)
//{
//  int ncells = data.ncells;
//  int nmats = data.nmats;
//  std::vector<double>& Volfrac = data.Volfrac;
//  std::vector<double>& Densityfrac = data.Densityfrac;
//  const std::vector<double>& cen = data.cen;
//  const std::vector<int>& nnbrs = data.nnbrs;
//  const std::vector<int>& nbrs = data.nbrs;
//
//
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
//void average_material_density_over_cell_nbr_mm(MultiMat& mm, Robey_data& data)
//{
//  cout << "-- Calculating avg material density over cell neighbor, Multimat version--" << endl;
//
//  const std::vector<double>& cen = data.cen;
//  const std::vector<int>& nnbrs = data.nnbrs;
//  const std::vector<int>& nbrs = data.nbrs;
//
//  int ncells = mm.getNumberOfCells();
//  int nmats = mm.getNumberOfMaterials();
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


int main(int argc, char** argv)
{
  //Robey_data data;

  ////Set-up the multimat class
  //MultiMat mm;
  //mm.setNumberOfMat(data.nmats);
  //mm.setNumberOfCell(data.ncells);
  //mm.setCellMatRel(data.Volfrac_bool);

  ////Setting field data in terms of slam
  //mm.addField<>("Densityfrac"    , FieldMapping::PER_CELL_MAT, &data.Densityfrac_sparse[0]);
  //mm.addField<>("Vol"            , FieldMapping::PER_CELL    , &data.Vol[0]);
  //mm.addField<>("Volfrac"        , FieldMapping::PER_CELL_MAT, &data.Volfrac_sparse[0]);
  //mm.addField<>("Temperaturefrac", FieldMapping::PER_CELL_MAT, &data.Temperaturefrac_sparse[0]);
  //mm.addField<>("Pressurefrac"   , FieldMapping::PER_CELL_MAT, &data.Pressurefrac_sparse[0]);
  //mm.addField<>("nmatconsts"     , FieldMapping::PER_MAT     , &data.nmatconsts[0]);

  ////printself and check
  //mm.printSelf();
  //printf("IsValid: %d\n\n", mm.isValid());

  ////Run the examples
  //average_density_cell_dom(data);
  //average_density_cell_dom_mm(mm);
  //average_density_mat_dom_mm(mm);

  //average_density_cell_dom_with_if(data);
  //average_density_cell_dom_with_if_mm(mm);

  //calculate_pressure(data);
  //calculate_pressure_mm(mm);

  //average_material_density_over_cell_nbr(data);
  //average_material_density_over_cell_nbr_mm(mm, data);

  test_code();

  cout << "\nPress Enter to terminate...";
  cin.ignore();
}


void test_code() {

  //Set-up
  bool use_sparse = false;

  MultiMat mm(DataLayout::CELL_CENTRIC, use_sparse ? SparcityLayout::SPARSE : SparcityLayout::DENSE);
  
  int nmats = 50;
  int ncells = 2000;
  int ncomp = 2;

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
    for (int s = 0; s < ncomp; ++s) {
      cell_arr1[i*ncomp + s] = (double)i * 2.0 + s * 0.01;
      c_sum += cell_arr1[i*ncomp + s];
    }
  }
    
  std::vector<double> cellmat_arr;
  cellmat_arr.resize( (use_sparse ? nfilled : nmats * ncells) *ncomp);
  double x_sum = 0;
  for (unsigned int i = 0; i < cellmat_arr.size()/ncomp; i++) {
    if (use_sparse || cellMatRel[i]) {
      for (int s = 0; s < ncomp; ++s) {
        cellmat_arr[i*ncomp + s] = (double)i * 1.1 + s * 0.001;
        x_sum += cellmat_arr[i*ncomp + s];
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
      for (int s = 0; s < map.stride(); s++)
      {
        double val = map(i, s);                               //<----
        assert(val == map[ i*map.numComp() + s] );    //1d bracket access
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
      MultiMat::SubField<double>& submap = map2d(i);
      assert(rel_set.size() == submap.size());
      for (int k = 0; k < submap.size(); k++)
      {
        int idx = submap.index(k);  //mat id
        int idx2 = rel_set[k]; //another way to get mat id
        assert(idx == idx2); 

        for (int c = 0; c < submap.numComp(); ++c) {
          double val = submap.value(k, c);                 //<----------
          assert( val == submap(k,c) );                    //operator () access
          assert( val == submap[ k*submap.stride()+c ] );  //1d bracket access
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
          for (int s = 0; s < map.numComp(); ++s) {
            double val = map[indexSet[j]*map.numComp() + s];   //<-----
            //assert(val == map(indexSet.at(j), s)); //if 1dMap is used, this is also valid
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
      for(int s=0; s<iter.numComp(); ++s)
        sum += iter(s);              //<----
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
      MultiMat::SubField<double>& submap = map2d(i);
      for (auto iter = submap.begin(); iter != submap.end(); iter++) 
      {
        for (int s = 0; s < map2d.stride(); ++s)
        {
          sum += iter(s);                  //<----
          assert(iter(s) == iter.value(s));  //another way to get the value
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
      for (MultiMat::SubField<double>::SubmapIterator iter = map2d.begin(i); iter != map2d.end(i); iter++)
      {
        for (int s = 0; s < map2d.numComp(); ++s)
        {
          sum += iter(s);          //<----
          assert(iter(s) == iter.value(s));  //another way to get the value
          int idx = iter.index();
        }
        assert(iter(0) == *iter); //2 ways to get the first component value
        assert(iter(0) == iter.value()); 
      }
    }
  }
  cout << end_timer() << "\t";
  assert(x_sum == sum);

  
  cout << endl;
}
