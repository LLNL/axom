
#include "multimat_example.hpp"

#include "helper.hpp"
#include <ctime>



using namespace std;
using namespace axom::multimat;

using MM_doubleArrType = MultiMatArray<double>;

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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    Average density with fractional densities
//		Copied from Robey's code for Cell-Dominant Full Matrix Data Structure line 257
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void average_density_cell_dom(int ncells, int nmats, std::vector<double>& Volfrac, std::vector<double>& Densityfrac, vector<double>& Vol) {
	cout << "-- Averaging Density cell-dominant array-access --" << endl;
	vector<double> Density_average(ncells);

	double time_sum = 0;
	for (int iter = 0; iter < ITERMAX; iter++) {
		std::clock_t start;
		double duration;
		start = std::clock();

		for (int ic = 0; ic < ncells; ic++) {
			double density_ave = 0.0;
			for (int m = 0; m < nmats; m++) {
				density_ave += Densityfrac[ic*nmats + m] * Volfrac[ic*nmats + m];
			}
			Density_average[ic] = density_ave / Vol[ic];
		}

		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		time_sum += duration;

		//Do a check for calculation correctness
		double sum = 0;
		for (auto da : Density_average)
			sum += da;
		cout << "iter" << iter << "   sum of values: " << sum << endl;
	}

	double act_perf = time_sum / (double)ITERMAX;
	printf("Average Density of mixed material cells    compute time is %lf secs\n", act_perf);

}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    Average density with fractional densities
//		Same as the function above, but modified to use MultiMat class
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void average_density_cell_dom_mm(MultiMat& mm) {
	cout << "-- Averaging Density cell-dominant mm-version --" << endl;
	int ncells = mm.m_ncells;
	int nmats = mm.m_nmats;
	MM_doubleArrType* Densityfrac = MM_CAST_TO(double, mm.getFieldArray("Densityfrac"));
	MM_doubleArrType* Volfrac = MM_CAST_TO(double, mm.getFieldArray("Volfrac"));
	MM_doubleArrType* Vol = MM_CAST_TO(double, mm.getFieldArray("Vol"));

	vector<double> Density_average(ncells);

	double time_sum = 0;

	for (int iter = 0; iter < ITERMAX; iter++) {
		std::clock_t start;
		double duration;
		start = std::clock();

		for (int ic = 0; ic < ncells; ic++) {
			double density_ave = 0.0;
			for (int m = 0; m < nmats; m++) {
				density_ave += Densityfrac->getValue(ic, m) *  Volfrac->getValue(ic, m);
			}
			Density_average[ic] = density_ave / Vol->getValue(ic);
		}

		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		time_sum += duration;

		//Do a check for calculation correctness
		double sum = 0;
		for (auto da : Density_average)
			sum += da;
		cout << "iter" << iter << "   sum of values: " << sum << endl;
	}
	double act_perf = time_sum / (double)ITERMAX;
	printf("Average Density of mixed material cells    compute time is %lf secs\n", act_perf);
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    Average density with fractional densities
//        Material-first loop
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void average_density_mat_dom_mm(MultiMat& mm) {
	cout << "-- Averaging Density mat-dominant mm-version --" << endl;
	int ncells = mm.m_ncells;
	int nmats = mm.m_nmats;
	MM_doubleArrType* Densityfrac = MM_CAST_TO(double, mm.getFieldArray("Densityfrac"));;
	MM_doubleArrType* Volfrac = MM_CAST_TO(double, mm.getFieldArray("Volfrac"));
	MM_doubleArrType* Vol = MM_CAST_TO(double, mm.getFieldArray("Vol"));

	double time_sum = 0;

	for (int iter = 0; iter < ITERMAX; iter++) {
		vector<double> Density_average(ncells);

		std::clock_t start;
		double duration;
		start = std::clock();

		for (int m = 0; m < nmats; m++) {
			for (int ic = 0; ic < ncells; ic++) {
				Density_average[ic] += Densityfrac->getValue(ic, m) * Volfrac->getValue(ic, m);
			}
		}
		for (int ic = 0; ic < ncells; ic++) {
			Density_average[ic] /= Vol->getValue(ic);
		}

		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		time_sum += duration;

		//Do a check for calculation correctness
		double sum = 0;
		for (auto da : Density_average)
			sum += da;
		cout << "iter" << iter << "   sum of values: " << sum << endl;
	}
	double act_perf = time_sum / (double)ITERMAX;
	printf("Average Density of mixed material cells    compute time is %lf secs\n", act_perf);
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    Average density with fractional densities
//		Same as the function above, but modified to use MultiMat class
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void average_density_cell_dom_iter_mm(MultiMat& mm) {
	cout << "-- Averaging Density cell-dominant using iterator mm-version --" << endl;
	//TODO when mm class has iterators
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    Average density with fractional densities with if test
//		Copied from Robey's code for Cell-Dominant Full Matrix Data Structure
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void average_density_cell_dom_with_if(int ncells, int nmats, std::vector<double>& Volfrac, std::vector<double>& Densityfrac, vector<double>& Vol)
{
	cout << "-- Averaging Density with if --" << endl;
	vector<double> Density_average(ncells, 0);
	double time_sum = 0;

	for (int iter = 0; iter < ITERMAX; iter++) {
		std::clock_t start;
		double duration;
		start = std::clock();

		for (int ic = 0; ic < ncells; ic++) {
			double density_ave = 0.0;
			for (int m = 0; m < nmats; m++) {
				if (Volfrac[ic*nmats + m] > 0.0) {
					density_ave += Densityfrac[ic*nmats + m] * Volfrac[ic*nmats + m];
				}
			}
			Density_average[ic] = density_ave / Vol[ic];
		}

		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		time_sum += duration;

		//Do a check for calculation correctness
		double sum = 0;
		for (auto da : Density_average)
			sum += da;
		cout << "iter" << iter << "   sum of values: " << sum << endl;
	}

	double act_perf = time_sum / (double)ITERMAX;
	printf("Average Density of frac with if            compute time is %lf secs\n", act_perf);

}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    Average density with fractional densities with if test
//          Same as the function above, but modified to use MultiMat class
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void average_density_cell_dom_with_if_mm(MultiMat& mm) {
	cout << "-- Averaging Density with if, mm-version --" << endl;
	int ncells = mm.m_ncells;
	int nmats = mm.m_nmats;
	MM_doubleArrType* Densityfrac = MM_CAST_TO(double, mm.getFieldArray("Densityfrac"));;
	MM_doubleArrType* Volfrac = MM_CAST_TO(double, mm.getFieldArray("Volfrac"));
	MM_doubleArrType* Vol = MM_CAST_TO(double, mm.getFieldArray("Vol"));

	vector<double> Density_average(ncells, 0);

	double time_sum = 0;

	for (int iter = 0; iter < ITERMAX; iter++) {
		std::clock_t start;
		double duration;
		start = std::clock();

		for (int ic = 0; ic < ncells; ic++) {
			double density_ave = 0.0;
			for (int m = 0; m < nmats; m++) {
				if (Volfrac->getValue(ic, m) > 0.0) {
					density_ave += Densityfrac->getValue(ic, m) *  Volfrac->getValue(ic, m);
				}
			}
			Density_average[ic] = density_ave / Vol->getValue(ic);
		}

		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		time_sum += duration;

		//Do a check for calculation correctness
		double sum = 0;
		for (auto da : Density_average)
			sum += da;
		cout << "iter" << iter << "   sum of values: " << sum << endl;
	}

	double act_perf = time_sum / (double)ITERMAX;
	printf("Average Density of frac with if            compute time is %lf secs\n", act_perf);

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   Calculate pressure using ideal gas law
//		Copied from Robey's code for Cell-Dominant Full Matrix Data Structure
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calculate_pressure(int ncells, int nmats, std::vector<double>& Volfrac, std::vector<double>& Densityfrac,
	vector<double>& Temperaturefrac, vector<double>& nmatconsts)
{
	cout << "-- Calculating pressure --" << endl;
	vector<double> Pressurefrac(ncells*nmats, 0);

	double time_sum = 0;

	for (int iter = 0; iter< ITERMAX; iter++) {
		std::clock_t start;
		double duration;
		start = std::clock();

		for (int ic = 0; ic < ncells; ic++) {
			for (int m = 0; m < nmats; m++) {
				if (Volfrac[ic*nmats + m] > 0.) {
					Pressurefrac[ic*nmats + m] = (nmatconsts[m] * Densityfrac[ic*nmats + m] * Temperaturefrac[ic*nmats + m]) / (Volfrac[ic*nmats + m]);
				}
				else {
					Pressurefrac[ic*nmats + m] = 0.0;
				}
			}
		}

		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		time_sum += duration;

		//Do a check for calculation correctness
		double sum = 0;
		for (auto da : Pressurefrac)
			sum += da;
		cout << "iter" << iter << "   sum of values: " << sum << endl;
	}
	double act_perf = time_sum / (double)ITERMAX;
	printf("Pressure Calculation of mixed material cells with if compute time is %lf msecs\n", act_perf);

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   Calculate pressure using ideal gas law
//          Same as the function above, but modified to use MultiMat class
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calculate_pressure_mm(MultiMat& mm)
{
	cout << "-- Calculating pressure, mm-version --" << endl;
	int ncells = mm.m_ncells;
	int nmats = mm.m_nmats;
	MM_doubleArrType* Densityfrac = MM_CAST_TO(double, mm.getFieldArray("Densityfrac"));;
	MM_doubleArrType* Volfrac = MM_CAST_TO(double, mm.getFieldArray("Volfrac"));
	//	MM_doubleArrType* Vol = MM_CAST_TO(double, mm.getFieldArray("Vol"));
	MM_doubleArrType* Temperaturefrac = MM_CAST_TO(double, mm.getFieldArray("Temperaturefrac"));

	vector<double> Pressurefrac(ncells*nmats, 0);

	double time_sum = 0;

	for (int iter = 0; iter< ITERMAX; iter++) {
		std::clock_t start;
		double duration;
		start = std::clock();

		for (int ic = 0; ic < ncells; ic++) {
			for (int m = 0; m < nmats; m++) {
				if (Volfrac->getValue(ic, m) > 0.) {
					Pressurefrac[ic*nmats + m] = ( /*nmatconsts[m]*/ 5.0 * Densityfrac->getValue(ic, m)
						* Temperaturefrac->getValue(ic, m)) / Volfrac->getValue(ic, m);
				}
				else {
					Pressurefrac[ic*nmats + m] = 0.0;
				}
			}
		}

		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		time_sum += duration;

		//Do a check for calculation correctness
		double sum = 0;
		for (auto da : Pressurefrac)
			sum += da;
		cout << "iter" << iter << "   sum of values: " << sum << endl;
	}
	double act_perf = time_sum / (double)ITERMAX;
	printf("Pressure Calculation of mixed material cells with if compute time is %lf msecs\n", act_perf);

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    Average material density over neighborhood of each cell
//		Copied from Robey's code for Cell-Dominant Full Matrix Data Structure
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void average_material_density_over_cell_nbr(int ncells, int nmats, std::vector<double>& Volfrac,
	std::vector<double>& Densityfrac, //vector<double>& Temperaturefrac, vector<double>& nmatconsts
	const std::vector<double>& cen, const std::vector<int>& nnbrs, const std::vector<int>& nbrs)
{
	cout << "-- Calculating avg material density over cell neighbor --" << endl;

	double time_sum = 0;

	for (int iter = 0; iter< ITERMAX; iter++) {
		vector<double> MatDensity_average(ncells*nmats, 0);

		std::clock_t start;
		double duration;
		start = std::clock();

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
						if (Volfrac[ic*nmats + m] > 0.0) {
							MatDensity_average[ic*nmats + m] += Densityfrac[ic*nmats + m] / dsqr[n];
							nnm++;
						}
					}
					MatDensity_average[ic*nmats + m] /= nnm;
				}
				else {
					MatDensity_average[ic*nmats + m] = 0.0;
				}
			}
		}

		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		time_sum += duration;

		//Do a check for calculation correctness
		double sum = 0;
		for (auto da : MatDensity_average)
			sum += da;
		cout << "iter" << iter << "   sum of values: " << sum << endl;
	}

	double act_perf = time_sum / (double)ITERMAX;

	printf("Average Material Density            compute time is %lf msecs\n", act_perf);

}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    Average material density over neighborhood of each cell
//		Same as the function above, but modified to use MultiMat class
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void average_material_density_over_cell_nbr_mm(MultiMat& mm, const std::vector<double>& cen, const std::vector<int>& nnbrs, const std::vector<int>& nbrs)
{
	cout << "-- Calculating avg material density over cell neighbor, Multimat version--" << endl;

	int ncells = mm.m_ncells;
	int nmats = mm.m_nmats;
	MM_doubleArrType* Densityfrac = MM_CAST_TO(double, mm.getFieldArray("Densityfrac"));;
	MM_doubleArrType* Volfrac = MM_CAST_TO(double, mm.getFieldArray("Volfrac"));
	//MM_doubleArrType* Vol = MM_CAST_TO(double, mm.getFieldArray("Vol"));
	//MM_doubleArrType* Temperaturefrac = MM_CAST_TO(double, mm.getFieldArray("Temperaturefrac"));

	double time_sum = 0;

	for (int iter = 0; iter< ITERMAX; iter++) {
		vector<double> MatDensity_average(ncells*nmats, 0);

		std::clock_t start;
		double duration;
		start = std::clock();

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
				if (Volfrac->getValue(ic, m)> 0.0) {
					int nnm = 0;         // number of nbrs with this material
					for (int n = 0; n < nn; n++) {
						int jc = cnbrs[n];
						if (Volfrac->getValue(ic, m) > 0.0) {
							MatDensity_average[ic*nmats + m] += Densityfrac->getValue(ic, m) / dsqr[n];
							nnm++;
						}
					}
					MatDensity_average[ic*nmats + m] /= nnm;
				}
				else {
					MatDensity_average[ic*nmats + m] = 0.0;
				}
			}
		}

		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		time_sum += duration;

		//Do a check for calculation correctness
		double sum = 0;
		for (auto da : MatDensity_average)
			sum += da;
		cout << "iter" << iter << "   sum of values: " << sum << endl;
	}

	double act_perf = time_sum / (double)ITERMAX;

	printf("Average Material Density            compute time is %lf msecs\n", act_perf);

}




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
	MM_doubleArrType* MMArr_densityfrac     = mm.newFieldArray<>("Densityfrac"    , PER_CELL_MAT, &Densityfrac_sparse[0]);
	MM_doubleArrType* MMArr_vol             = mm.newFieldArray<>("Vol"            , PER_CELL    , &Vol[0]);
	MM_doubleArrType* MMArr_volfrac         = mm.newFieldArray<>("Volfrac"        , PER_CELL_MAT, &Volfrac_sparse[0]);
	MM_doubleArrType* MMArr_temperaturefrac = mm.newFieldArray<>("Temperaturefrac", PER_CELL_MAT, &Temperaturefrac_sparse[0]);
	MM_doubleArrType* MMArr_pressurefrac    = mm.newFieldArray<>("Pressurefrac"   , PER_CELL_MAT, &Pressurefrac[0]);
	MM_doubleArrType* MMArr_nmatconsts      = mm.newFieldArray<>("nmatconsts"     , PER_MAT     , &nmatconsts[0]);

	//printself and check
	mm.printSelf();
	printf("IsValid: %d\n\n", mm.isValid());

	//data value check, to make sure the values in multimat class is set-up correctly
	for (int ci = 0; ci<mm.m_ncells; ci++) {
		for (int mi = 0; mi < mm.m_nmats; mi++) {
			double v = MMArr_densityfrac->getValue(ci, mi);
			assert(v == Densityfrac[ci*mm.m_nmats + mi]);

			v = MMArr_volfrac->getValue(ci, mi);
			assert(v == Volfrac[ci*mm.m_nmats + mi]);

			v = MMArr_temperaturefrac->getValue(ci, mi);
			assert(v == Temperaturefrac[ci*mm.m_nmats + mi]);

			v = MMArr_pressurefrac->getValue(ci, mi);
			assert(v == Pressurefrac[ci*mm.m_nmats + mi]);
		}
		assert(MMArr_vol->getValue(ci) == Vol[ci]);
	}
	for (int mi = 0; mi < mm.m_nmats; mi++) {
		double v = MMArr_nmatconsts->getValue(mi);
		assert(v == nmatconsts[mi]);
	}


	//Run the examples
	average_density_cell_dom(ncells, nmats, Volfrac, Densityfrac, Vol);
	average_density_cell_dom_mm(mm);
	average_density_mat_dom_mm(mm);

	average_density_cell_dom_with_if(ncells, nmats, Volfrac, Densityfrac, Vol);
	average_density_cell_dom_with_if_mm(mm);

	calculate_pressure(ncells, nmats, Volfrac, Densityfrac, Temperaturefrac, nmatconsts);
	calculate_pressure_mm(mm);

	average_material_density_over_cell_nbr(ncells, nmats, Volfrac, Densityfrac, cen, nnbrs, nbrs);
	average_material_density_over_cell_nbr_mm(mm, cen, nnbrs, nbrs);

	test_code();

	cout << "\nPress Enter to terminate...";
	cin.ignore();
}


void test_code() {
	//Set-up
	MultiMat mm;
	int nmats = 50;
	int ncells = 2000;
	mm.setNumberOfMat(nmats);
	mm.setNumberOfCell(ncells);
	std::vector<bool> cellMatRel(nmats * ncells, true);
	mm.setCellMatRel(cellMatRel);

	//add arrays
	std::vector<double> cell_arr1(ncells);
	double c_sum = 0;
	for (int i = 0; i < cell_arr1.size(); i++) {
		cell_arr1[i] = (double)i * 2.0;
		c_sum += cell_arr1[i];
	}
		
	auto MMarr_cell = mm.newFieldArray<>("Cell Array", PER_CELL, &cell_arr1[0]);
	MMarr_cell->setValue(0, 123);
	c_sum += 123;

	std::vector<double> cellmat_arr(nmats * ncells);
	double x_sum = 0;
	for (int i = 0; i < cellmat_arr.size(); i++) {
		cellmat_arr[i] = (double)i * 1.1;
		x_sum += cellmat_arr[i];
	}
	auto MMarr_cellmat = mm.newFieldArray<>("CellMat Array", PER_CELL_MAT, &cellmat_arr[0]);

	//Different accessing methods ...

	// ------------ templated get -----------
	printf("\nTemplate get function\n");
	double sum = 0;
	start_timer();
	for (int i = 0; i < mm.m_ncells; i++) {
		sum += mm.getFieldValue<double>("Cell Array", i); //<----
	}
	assert(c_sum == sum);
	cout << end_timer() << " ";

	sum = 0;
	start_timer();
	for (int i = 0; i < mm.m_ncells; i++) {
		for (int m = 0; m < mm.m_nmats; m++) {
			sum += mm.getFieldValue<double>("CellMat Array", i, m); //<----
		}
	}
	cout << end_timer() << " ";
	assert(x_sum == sum);


	// ------- Accessing using MMArray ----------
	printf("\nAccess through MMArray object\n");
	sum = 0;
	start_timer();
	auto arr3 = mm.getFieldArray("Cell Array");
	MultiMatArray<double>* arr4 = dynamic_cast<MultiMatArray<double>*>(arr3); //dynamic cast
	for (int i = 0; i < mm.m_ncells; i++) {
		sum += arr4->getValue(i);                  //<----
	}
	cout << end_timer() << " ";
	assert(c_sum == sum);

	sum = 0;
	MultiMatArray<double>* arr5 = dynamic_cast<MultiMatArray<double>*>(mm.getFieldArray("CellMat Array")); //dynamic cast
	start_timer();
	for (int i = 0; i < mm.m_ncells; i++) {
		//int matcount = arr5->getNumMatInCell(i);
		for (int m = 0; m < mm.m_nmats /*matcount*/; m++) {
			sum += arr5->getValue(i, m);          //<---- warning: extra for-loop in the call!
		}
	}
	cout << end_timer() << " ";
	assert(x_sum == sum);


	//Can also use macro to convert to typed MMARray
	MultiMatArray<double>* arr6 = MM_CAST_TO(double, mm.getFieldArray("Cell Array"));
	assert(arr4 == arr6);
	MultiMatArray<double>* arr7 = MM_CAST_TO(double, mm.getFieldArray("CellMat Array"));
	assert(arr5 == arr7);

	// --------- returning SLAM map -----------
	printf("\nWith SLAM map\n");
	sum = 0;
	start_timer();
	const MultiMatArray<double>::MapType Map = MMarr_cell->getMap(); //TODO reference
	for (int i = 0; i < mm.m_ncells; i++) {
		sum += Map[i];                               //<----
	}
	cout << end_timer() << " ";
	assert(c_sum == sum);

	/* //currently don't have a way to access cell-mat relation map through map type
	sum = 0;
	const MultiMatArray<double>::MapType Map2 = arr2->getMap();
	start_timer();
	for (int i = 0; i < mm.m_ncells; i++) {
		for (int m = 0; m < mm.m_nmats; m++) {
			int j = set.get(i,m);
			sum += Map[j];
		}
	}
	cout << end_timer() << " ";
	assert(x_sum == sum);
	*/

	// ---------- using iterator -------------
	printf("\nWith itereators\n");
	sum = 0;
	start_timer();
	for (MultiMatArray<double>::iterator a = MMarr_cell->begin(); a != MMarr_cell->end(); a++) {
		sum += *a;              //<----
	}
	cout << end_timer() << " ";
	assert(c_sum == sum);

	sum = 0;
	start_timer();
	for (MultiMatArray<double>::iterator a = MMarr_cellmat->begin(); a != MMarr_cellmat->end(); a++) {
		sum += *a;            //<----
	}
	cout << end_timer() << " ";
	assert(x_sum == sum);


	cout << endl;
}