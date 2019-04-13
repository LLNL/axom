// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * Set-up the test multi-material data.
 * From Robey's code for cell-dominant full matrix
 */

/*
struct full 
{
  vector<double> Volfrac;
  float filled_percentage;

  vector<double> Vol;
  vector<double> Density;
  vector<double> Temperature;
  vector<double> Pressure;
  vector<double> Densityfrac;
  vector<double> Temperaturefrac;
  vector<double> Pressurefrac;
};


struct compact
{
  vector<double> Volfrac_fullcc;
  float filled_percentage;

    int*& imaterial
    int*& nmaterials
    double*& Vol
    double*& Density
    double*& Temperature
    double*& Pressure
    int*& imaterialfrac
    int*& nextfrac,
    int*& frac2cell, 
    double*& Volfrac
    double*& Densityfrac
    double*& Temperaturefrac
    double*& Pressurefrac
};
*/


void make_other_field_data_celldom(int ncells, int nmats,
                           std::vector<double>& i_Volfrac_CD,
                           std::vector<double>& o_Volfrac,
                           std::vector<double>& o_Vol,
                           std::vector<double>& o_Densityfrac,
                           std::vector<double>& o_Temperaturefrac,
                           std::vector<double>& o_Pressurefrac,
                           std::vector<double>& o_Volfrac_sparse,
                           std::vector<double>& o_Densityfrac_sparse,
                           std::vector<double>& o_Temperaturefrac_sparse,
                           std::vector<double>& o_Pressurefrac_sparse,
                           std::vector<int>& o_begin_idx,
                           std::vector<int>& o_col_idx
  )
{
  //Set up dense data
  o_Vol.resize(ncells, 0);
  o_Volfrac.resize(ncells*nmats, 0);
  o_Densityfrac.resize(ncells*nmats, 0);
  o_Temperaturefrac.resize(ncells*nmats, 0);
  o_Pressurefrac.resize(ncells*nmats, 0);
  int nnz = 0;
  for (int ic = 0 ; ic < ncells ; ic++)
  {
    o_Vol[ic] = 0.0;
    for (int m = 0 ; m < nmats ; m++)
    {
      o_Volfrac[ic*nmats + m] = i_Volfrac_CD[ic*nmats + m];
      if (i_Volfrac_CD[ic*nmats + m] > 0.0)
      {
        nnz += 1;
        o_Densityfrac[ic*nmats + m] = 2.0;
        o_Temperaturefrac[ic*nmats + m] = 0.5;
      }
      else
      {
        o_Densityfrac[ic * nmats + m] = 0.0;
        o_Temperaturefrac[ic * nmats + m] = 0.0;
      }
      o_Pressurefrac[ic * nmats + m] = 0.0;

      o_Vol[ic] += i_Volfrac_CD[ic*nmats + m];
    }
  }

  //fill in sparse data
  o_begin_idx.resize(ncells + 1);
  o_col_idx.resize(nnz);
  int ii = 0;
  for (int ic = 0; ic < ncells; ic++)
  {
    o_begin_idx[ic] = ii;
    for (int m = 0; m < nmats; m++)
    {
      if (i_Volfrac_CD[ic * nmats + m] > 0.0)
      {
        o_col_idx[ii] = m;
        ii++;
      }
    }
  }
  o_begin_idx[ncells] = ii;
  assert(ii == nnz);
  
  o_Densityfrac_sparse.resize(nnz);
  o_Volfrac_sparse.resize(nnz);
  o_Temperaturefrac_sparse.resize(nnz);
  o_Pressurefrac_sparse.resize(nnz);

  ii = 0;
  for (int ci = 0; ci < ncells; ci++)
  {
    for (int mi = 0; mi < nmats; mi++)
    {
      if (o_Volfrac[ci * nmats + mi] != 0)
      {
        o_Densityfrac_sparse[ii] = o_Densityfrac[ci * nmats + mi];
        o_Volfrac_sparse[ii] = o_Volfrac[ci * nmats + mi];
        o_Temperaturefrac_sparse[ii] = o_Temperaturefrac[ci * nmats + mi];
        o_Pressurefrac_sparse[ii] = o_Pressurefrac[ci * nmats + mi];
        ii++;
      }
    }
  }
  assert(ii == nnz);
}




void make_other_field_data_matdom(int ncells, int nmats,
                          std::vector<double>& i_Volfrac_CD,
                          std::vector<double>& o_Volfrac,
                          std::vector<double>& o_Vol,
                          std::vector<double>& o_Densityfrac,
                          std::vector<double>& o_Temperaturefrac,
                          std::vector<double>& o_Pressurefrac,
                          std::vector<double>& o_Volfrac_sparse,
                          std::vector<double>& o_Densityfrac_sparse,
                          std::vector<double>& o_Temperaturefrac_sparse,
                          std::vector<double>& o_Pressurefrac_sparse,
                          std::vector<int>& o_begin_idx,
                          std::vector<int>& o_col_idx
)
{
  //Set up dense data
  o_Vol.resize(ncells, 0);
  o_Volfrac.resize(ncells * nmats, 0);
  o_Densityfrac.resize(ncells * nmats, 0);
  o_Temperaturefrac.resize(ncells * nmats, 0);
  o_Pressurefrac.resize(ncells * nmats, 0);
  int nnz = 0;
  for (int ic = 0; ic < ncells; ic++)
  {
    o_Vol[ic] = 0.0;
  }

  for (int m = 0; m < nmats; m++)
  {
    for (int ic = 0; ic < ncells; ic++)
    {
      o_Volfrac[m*ncells + ic] = i_Volfrac_CD[ic * nmats + m];
      if (i_Volfrac_CD[ic * nmats + m] > 0.0)
      {
        nnz += 1;
        o_Densityfrac[m*ncells + ic] = 2.0;
        o_Temperaturefrac[m * ncells + ic] = 0.5;
      }
      else 
      {
        o_Densityfrac[m * ncells + ic] = 0;
        o_Temperaturefrac[m * ncells + ic] = 0;
      }
      o_Pressurefrac[m * ncells + ic] = 0.0;
      o_Vol[ic] += i_Volfrac_CD[ic * nmats + m];
    }
  }

  //fill in sparse data, CSR layout
  o_begin_idx.resize(nmats + 1);
  o_col_idx.resize(nnz);
  int ii = 0;
  for (int m = 0; m < nmats; m++)
  {
    o_begin_idx[m] = ii;
    for (int ic = 0; ic < ncells; ic++)
    {
      if (i_Volfrac_CD[ic * nmats + m] > 0.0)
      {
        o_col_idx[ii] = ic;
        ii++;
      }
    }
  }
  o_begin_idx[nmats] = ii;
  assert(ii == nnz);

  o_Densityfrac_sparse.resize(nnz);
  o_Volfrac_sparse.resize(nnz);
  o_Temperaturefrac_sparse.resize(nnz);
  o_Pressurefrac_sparse.resize(nnz);

  ii = 0;
  for (int mi = 0; mi < nmats; mi++)
  {
    for (int ci = 0; ci < ncells; ci++)
    {
      if (o_Volfrac[mi * ncells + ci] != 0)
      {
        o_Densityfrac_sparse[ii] = o_Densityfrac[mi * ncells + ci];
        o_Volfrac_sparse[ii] = o_Volfrac[mi * ncells + ci];
        o_Temperaturefrac_sparse[ii] = o_Temperaturefrac[mi * ncells + ci];
        o_Pressurefrac_sparse[ii] = o_Pressurefrac[mi * ncells + ci];
        ii++;
      }
    }
  }
  assert(ii == nnz);
}



void read_vol_frac_matrix_file(std::string filename, int& ncells, int& nmats,
                               std::vector<double> &Volfrac,
                               float& filled_percentage)
{
  using namespace std;

  ncells = 1000000;

  int status;
  FILE* fp;
  fp = fopen(filename.c_str(), "r");
  if (!fp)
  {
    fprintf(stderr, "unable to read volume fractions from file \"%s\"\n",
            filename.c_str());
    exit(-1);
  }

  status = fscanf(fp, "%d", &nmats);
  if (status < 0)
  {
    printf("error in read at line %d\n", __LINE__);
    exit(1);
  }

  Volfrac.resize(ncells * nmats, 0.0);

  char matname[256];
  for (int m = 0 ; m < nmats ; m++)
  {
    status = fscanf(fp, "%s", matname);            // read and discard
    if (status < 0)
    {
      printf("error in read at line %d\n", __LINE__);
      exit(1);
    }
  }

  double VolTotal = 0.0;
  int filled_count = 0;
  int mixed_cell_count = 0;
  int mixed_frac_count = 0;
  int pure_frac_count = 0;
  int pure_cell_count = 0;
  int onematcell = 0;
  int twomatcell = 0;
  int threematcell = 0;
  int fourmatcell = 0;
  int fiveplusmatcell = 0;
  for (int ic = 0 ; ic < ncells ; ic++)
  {
    int mat_count = 0;
    for (int m = 0 ; m < nmats ; m++)
    {
      status = fscanf(fp, "%lf", &(Volfrac[ic*nmats + m]));
      if (status < 0)
      {
        printf("error in read at line %d\n", __LINE__);
        exit(1);
      }
      if (Volfrac[ic*nmats + m] > 0.0)
      {
        filled_count++;
        mat_count++;
      }
      VolTotal += Volfrac[ic*nmats + m];
    }
    if (mat_count >= 2)
    {
      mixed_cell_count++;
      mixed_frac_count += mat_count;
    }
    else
    {
      pure_frac_count++;
    }
    if (mat_count == 1)
      pure_cell_count++;
    if (mat_count == 1)
      onematcell++;
    if (mat_count == 2)
      twomatcell++;
    if (mat_count == 3)
      threematcell++;
    if (mat_count == 4)
      fourmatcell++;
    if (mat_count >= 5)
      fiveplusmatcell++;
  }
  fclose(fp);

  printf("Ratios to Full Data Structure\n");
  filled_percentage = (float)filled_count*100.0 / (float)(ncells*nmats);
  float sparsity_percentage = (float)(ncells*nmats - filled_count)*100.0 /
                              (float)(ncells*nmats);
  printf("Sparsity %lf percent/Filled %lf percent\n\n",
         sparsity_percentage, filled_percentage);

  printf("Ratios to Number of Cells\n");
  float pure_cell_percentage = (float)pure_cell_count*100.0 / (float)ncells;
  float mixed_cell_percentage = (float)mixed_cell_count*100.0 / (float)ncells;
  printf("Pure cell %lf percentage/Mixed material %lf percentage\n\n",
         pure_cell_percentage, mixed_cell_percentage);

  printf("Ratios to Mixed Material Data Structure\n");
  float mixed_material_sparsity_percentage = (float)mixed_frac_count*100.0 /
                                             (float)(mixed_cell_count*nmats);
  float mixed_material_filled_percentage =
    (float)(mixed_cell_count*nmats - mixed_frac_count)*100.0 /
    (float)(mixed_cell_count*nmats);
  printf(
    "Mixed material Sparsity %lf percent/Mixed material Filled %lf percent\n\n",
    mixed_material_sparsity_percentage, mixed_material_filled_percentage);

  printf("Vol Total %lf\n", VolTotal);
  printf("%f percent of the cells are filled\n",
         (float)filled_count*100.0 / (float)(ncells*nmats));
  printf("%f percent of the cells are mixed\n",
         (float)mixed_cell_count*100.0 / (float)ncells);
  printf("%f percent of the total are mixed\n",
         (float)mixed_frac_count*100.0 / (float)(ncells*nmats));
  printf("%f percent of the frac are mixed\n",
         (float)mixed_frac_count*100.0 / (float)(mixed_cell_count*nmats));
  printf("%f percent sparsity\n",
         (float)(ncells*nmats - mixed_frac_count)*100.0 /
         (float)(ncells*nmats));
  printf("%f percent of the frac are pure\n",
         (float)pure_frac_count*100.0 / (float)ncells);
  printf("1 matcell %d 2 matcell %d 3 matcell %d 4 matcell %d 5 matcell %d\n\n",
         onematcell, twomatcell, threematcell, fourmatcell, fiveplusmatcell);
  printf("Total cells %d\n\n",
         onematcell + 2 * twomatcell + 3 * threematcell + 4 * fourmatcell + 5 *
         fiveplusmatcell);

}

void get_vol_frac_matrix_rand(int& ncells, int& nmats,
                              std::vector<double> &Volfrac,
                              float& filled_percentage,
                              int default_ncell = 1000000,
                              int default_nmats = 50) {
//void get_vol_frac_matrix_rand(double **&Volfrac, float& filled_percentage) {
  using namespace std;

  ncells = default_ncell;
  nmats = default_nmats;

  Volfrac.resize(ncells * nmats, 0.0);

  vector<int> mf_rand(ncells, 0);
  srand(0);
  for (int ic = 0 ; ic < ncells ; ic++)
  {
    mf_rand[ic] =
      (int)((float)rand()*1000.0 / (float)((long long)RAND_MAX + 1));
  }

  double VolTotal = 0.0;
  int filled_count = 0;
  int mixed_cell_count = 0;
  int mixed_frac_count = 0;
  int pure_frac_count = 0;
  int pure_cell_count = 0;
  int onematcell = 0;
  int twomatcell = 0;
  int threematcell = 0;
  int fourmatcell = 0;
  int fiveplusmatcell = 0;
  for (int ic = 0 ; ic < ncells ; ic++)
  {
    int m1 =
      (int)((float)rand()*(float)nmats / (float)((long long)RAND_MAX + 1));
    if (m1 > nmats - 1)
      m1 = nmats - 1;
    Volfrac[ic*nmats+m1] = 1.0;
    int mf = mf_rand[ic];
    if (mf < 25)
    {
      int m2 =
        (int)((float)rand()*(float)nmats / (float)((long long)RAND_MAX + 1));
      while (m2 == m1)
      {
        m2 =
          (int)((float)rand()*(float)nmats / (float)((long long)RAND_MAX + 1));
      }
      int m3 =
        (int)((float)rand()*(float)nmats / (float)((long long)RAND_MAX + 1));
      while (m3 == m2 || m3 == m1)
      {
        m3 =
          (int)((float)rand()*(float)nmats / (float)((long long)RAND_MAX + 1));
      }
      int m4 =
        (int)((float)rand()*(float)nmats / (float)((long long)RAND_MAX + 1));
      while (m4 == m3 || m4 == m2 || m4 == m1)
      {
        m4 =
          (int)((float)rand()*(float)nmats / (float)((long long)RAND_MAX + 1));
      }
      if (m2 > nmats - 1)
        m2 = nmats - 1;
      if (m3 > nmats - 1)
        m3 = nmats - 1;
      if (m4 > nmats - 1)
        m4 = nmats - 1;
      Volfrac[ic*nmats + m1] = 0.4;
      Volfrac[ic*nmats + m2] = 0.3;
      Volfrac[ic*nmats + m3] = 0.2;
      Volfrac[ic*nmats + m4] = 0.1;
    }
    else if (mf < 75)
    {
      int m2 = (float)rand()*(float)nmats / (float)((long long)RAND_MAX + 1);
      while (m2 == m1)
      {
        m2 = (float)rand()*(float)nmats / (float)((long long)RAND_MAX + 1);
      }
      int m3 = (float)rand()*(float)nmats / (float)((long long)RAND_MAX + 1);
      while (m3 == m2 || m3 == m1)
      {
        m3 = (float)rand()*(float)nmats / (float)((long long)RAND_MAX + 1);
      }
      if (m2 > nmats - 1)
        m2 = nmats - 1;
      if (m3 > nmats - 1)
        m3 = nmats - 1;
      Volfrac[ic*nmats + m1] = 0.5;
      Volfrac[ic*nmats + m2] = 0.3;
      Volfrac[ic*nmats + m3] = 0.2;
    }
    else if (mf < 200)
    {
      int m2 = (float)rand()*(float)nmats / (float)((long long)RAND_MAX + 1);
      while (m2 == m1)
      {
        m2 = (float)rand()*(float)nmats / (float)((long long)RAND_MAX + 1);
      }
      if (m2 > nmats - 1)
        m2 = nmats - 1;
      Volfrac[ic*nmats + m1] = 0.5;
      Volfrac[ic*nmats + m2] = 0.5;
    }
    int mat_count = 0;
    for (int m = 0 ; m < nmats ; m++)
    {
      if (Volfrac[ic*nmats + m] > 0.0)
      {
        filled_count++;
        mat_count++;
      }
      VolTotal += Volfrac[ic*nmats + m];
    }
    if (mat_count >= 2)
    {
      mixed_cell_count++;
      mixed_frac_count += mat_count;
    }
    else
    {
      pure_frac_count++;
    }
    if (mat_count == 1)
      pure_cell_count++;
    if (mat_count == 1)
      onematcell++;
    if (mat_count == 2)
      twomatcell++;
    if (mat_count == 3)
      threematcell++;
    if (mat_count == 4)
      fourmatcell++;
    if (mat_count >= 5)
      fiveplusmatcell++;
  }

  printf("Ratios to Full Data Structure\n");
  filled_percentage = (float)filled_count*100.0 / (float)(ncells*nmats);
  float sparsity_percentage = (float)(ncells*nmats - filled_count)*100.0 /
                              (float)(ncells*nmats);
  printf("Sparsity %lf percent/Filled %lf percent\n\n",
         sparsity_percentage, filled_percentage);

  printf("Ratios to Number of Cells\n");
  float pure_cell_percentage = (float)pure_cell_count*100.0 / (float)ncells;
  float mixed_cell_percentage = (float)mixed_cell_count*100.0 / (float)ncells;
  printf("Pure cell %lf percentage/Mixed material %lf percentage\n\n",
         pure_cell_percentage, mixed_cell_percentage);

  printf("Ratios to Mixed Material Data Structure\n");
  float mixed_material_sparsity_percentage = (float)mixed_frac_count*100.0 /
                                             (float)(mixed_cell_count*nmats);
  float mixed_material_filled_percentage =
    (float)(mixed_cell_count*nmats - mixed_frac_count)*100.0 /
    (float)(mixed_cell_count*nmats);
  printf(
    "Mixed material Sparsity %lf percent/Mixed material Filled %lf percent\n\n",
    mixed_material_sparsity_percentage, mixed_material_filled_percentage);

  printf("Vol Total %lf\n",VolTotal);
  printf("%f percent of the cells are filled\n",
         (float)filled_count*100.0/(float)(ncells*nmats));
  printf("%f percent of the cells are mixed\n",
         (float)mixed_cell_count*100.0/(float)ncells);
  printf("%f percent of the total are mixed\n",
         (float)mixed_frac_count*100.0/(float)(ncells*nmats));
  printf("%f percent of the frac are mixed\n",
         (float)mixed_frac_count*100.0/(float)(mixed_cell_count*nmats));
  printf("%f percent sparsity\n",
         (float)(ncells*nmats-mixed_frac_count)*100.0/(float)(ncells*nmats));
  printf("%f percent of the frac are pure\n",
         (float)pure_frac_count*100.0/(float)ncells);
  printf("1 matcell %d 2 matcell %d 3 matcell %d 4 matcell %d 5 matcell %d\n\n",
         onematcell, twomatcell, threematcell, fourmatcell, fiveplusmatcell);
  printf("Total cells %d\n\n",
         onematcell+2*twomatcell+3*threematcell+4*fourmatcell+5*
         fiveplusmatcell);

}

void get_neighbors(int ncells, int nnbrs_max, std::vector<int>& num_nbrs,
                   std::vector<int>& nbrs) {
  int ncells1 = (int)sqrt(ncells);  // assumes ncells is a perfect square
  if (ncells1*ncells1 != ncells)
  {
    fprintf(stderr, "Number of cells in mesh is not a perfect square");
    exit(-1);
  }
  if (nnbrs_max != 8)
  {
    std::cout << "This neighbor generating code assumes 2D 8-neighbors. "
              << "You need to update this." << std::endl;
    exit(-1);
  }

  for (int i = 0 ; i < ncells1 ; i++)
  {
    for (int j = 0 ; j < ncells1 ; j++)
    {
      int c = i*ncells1 + j;
      int ilo = i == 0 ? i : i - 1;
      int jlo = j == 0 ? j : j - 1;
      int ihi = i == ncells1 - 1 ? i : i + 1;
      int jhi = j == ncells1 - 1 ? j : j + 1;
      int n = 0;
      for (int i1 = ilo ; i1 <= ihi ; i1++)
        for (int j1 = jlo ; j1 <= jhi ; j1++)
        {
          int c2 = i1*ncells1 + j1;
          if (c2 != c)
          {
            nbrs[c*nnbrs_max+n] = i1*ncells1 + j1;
            n++;
          }
        }
      num_nbrs[c] = n;
    }
  }

}

void get_centroids(int ncells, std::vector<double>& cen) {
  int ncells1 = (int)sqrt(ncells);  // assumes ncells is a perfect square
  if (ncells1*ncells1 != ncells)
  {
    fprintf(stderr, "Number of cells in mesh is not a perfect square");
    exit(-1);
  }

  // Assume domain is a unit square

  double XLO = 0.0, YLO = 0.0, XHI = 1.0, YHI = 1.0;
  double dx = (XHI - XLO) / ncells1, dy = (YHI - YLO) / ncells1;

  for (int i = 0 ; i < ncells1 ; i++)
  {
    for (int j = 0 ; j < ncells1 ; j++)
    {
      int c = i*ncells1 + j;
      cen[c*2] = XLO + i*dx;
      cen[c*2+1] = YLO + j*dy;
    }
  }
}
