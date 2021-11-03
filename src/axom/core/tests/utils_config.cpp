// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"

#ifdef AXOM_USE_OPENMP
  #include <omp.h>
#endif

#ifdef AXOM_USE_MFEM
  #include "mfem.hpp"
#endif

#include <algorithm>  // copy
#include <iostream>   // for cout
#include <sstream>    // for stringstream
#include <string>     // for C++ string
#include <vector>     // for STL vector

TEST(utils_config, axom_version)
{
  const int AXOM_MAJOR = AXOM_VERSION_MAJOR;
  const int AXOM_MINOR = AXOM_VERSION_MINOR;
  const int AXOM_PATCH = AXOM_VERSION_PATCH;
  EXPECT_TRUE(AXOM_MAJOR >= 0);
  EXPECT_TRUE(AXOM_MINOR >= 0);
  EXPECT_TRUE(AXOM_PATCH >= 0);

  const std::string AXOM_FULL = AXOM_VERSION_FULL;
  EXPECT_FALSE(AXOM_FULL.empty());

  std::cout << "AXOM_VERSION_FULL: " << AXOM_VERSION_FULL << std::endl;
  std::cout << "AXOM_VERSION_MAJOR: " << AXOM_VERSION_MAJOR << std::endl;
  std::cout << "AXOM_VERSION_MINOR: " << AXOM_VERSION_MINOR << std::endl;
  std::cout << "AXOM_VERSION_PATCH: " << AXOM_VERSION_PATCH << std::endl;

#ifdef AXOM_VERSION_EXTRA
  const std::string AXOM_EXTRA = AXOM_VERSION_EXTRA;
  EXPECT_FALSE(AXOM_EXTRA.empty());

  std::cout << "AXOM_VERSION_EXTRA: " << AXOM_VERSION_EXTRA << std::endl;
#endif
}

TEST(utils_config, config_libraries)
{
  // This test checks which libraries are available in the configuration

  std::cout << "Available libraries: " << std::endl;

  std::vector<std::string> libs;

#ifdef AXOM_USE_C2C
  libs.push_back("c2c");
#endif

#ifdef AXOM_USE_CLI11
  libs.push_back("CLI11");
#endif

#ifdef AXOM_USE_CONDUIT
  libs.push_back("conduit");
#endif

  libs.push_back(AXOM_CXX_STD);

#ifdef AXOM_USE_FMT
  libs.push_back("fmt");
#endif

#ifdef AXOM_USE_HDF5
  libs.push_back("hdf5");
#endif

#ifdef AXOM_USE_OPENMP
  libs.push_back("openmp");
#endif

#ifdef AXOM_USE_MFEM
  libs.push_back("mfem");
#endif

#ifdef AXOM_USE_MPI
  libs.push_back("mpi");
#endif

#ifdef AXOM_USE_MPIF_HEADERS
  libs.push_back("mpif headers");
#endif

#ifdef AXOM_USE_SPARSEHASH
  libs.push_back("sparsehash");
#endif

  libs.push_back("std::unordered_map");

  std::stringstream sstr;
  std::copy(libs.begin(),
            libs.end(),
            std::ostream_iterator<std::string>(sstr, "; "));
  std::cout << "\t{ " << sstr.str() << "}" << std::endl;

  EXPECT_TRUE(true);
}

TEST(utils_config, config_components)
{
  // This test checks which toolkit components are available in the
  // configuration

  std::cout << "Available components: " << std::endl;

  std::vector<std::string> comps;

  comps.push_back("core");

#ifdef AXOM_USE_MINT
  comps.push_back("mint");
#endif

#ifdef AXOM_USE_LUMBERJACK
  comps.push_back("lumberjack");
#endif

#ifdef AXOM_USE_PRIMAL
  comps.push_back("primal");
#endif

#ifdef AXOM_USE_QUEST
  comps.push_back("quest");
#endif

#ifdef AXOM_USE_SIDRE
  comps.push_back("sidre");
#endif

#ifdef AXOM_USE_SLAM
  comps.push_back("slam");
#endif

#ifdef AXOM_USE_SLIC
  comps.push_back("slic");
#endif

  std::stringstream sstr;
  std::copy(comps.begin(),
            comps.end(),
            std::ostream_iterator<std::string>(sstr, "; "));
  std::cout << "\t{ " << sstr.str() << "}" << std::endl;

  EXPECT_TRUE(true);
}

TEST(utils_config, config_openmp)
{
  // This test checks that the per-target OpenMP guards
  // in our configuration file 'axom/config.hpp' are working properly

#ifdef AXOM_USE_OPENMP
  std::cout << "OpenMP is available in this configuration." << std::endl;

  #pragma omp parallel
  {
    int thId = omp_get_thread_num();
    int thNum = omp_get_num_threads();
    int thMax = omp_get_max_threads();

  #pragma omp critical
    std::cout << "\tMy thread id is: " << thId << "\tNum threads is: " << thNum
              << "\tMax threads is: " << thMax << std::endl;
  }

#else
  std::cout << "OpenMP not available in this configuration." << std::endl;
#endif

  // Tests a simple reduction over integers
  // Note: Non-openmp configurations are compiled with a flag
  //       to ignore usage of unknown openmp pragmas
  const int N = 100;
  int sum = 0;
#pragma omp parallel for reduction(+ : sum)
  for(int i = 1; i <= N; ++i)
  {
    sum += i;
  }
  EXPECT_EQ(N * (N + 1) / 2, sum)
    << "Bad reduction on first " << N << " integers";
  std::cout << "Sum of first " << N << " numbers is: " << sum << std::endl;
}

#ifdef AXOM_USE_MFEM
TEST(utils_config, mfem_configuration)
{
  #ifdef MFEM_VERSION
  std::cout << "Using mfem version " << MFEM_VERSION_MAJOR
            << "."                        // major version
            << MFEM_VERSION_MINOR << "."  // minor version
            << MFEM_VERSION_PATCH         // patch level
            << std::endl;
  #endif  // MFEM_VERSION

  // Verify that this copy of mfem is configured appropriately with respect to MPI
  {
    bool axomHasMPI = false;
  #ifdef AXOM_USE_MPI
    axomHasMPI = true;
  #endif
    std::cout << "Axom is built " << (axomHasMPI ? "with" : "without") << " MPI"
              << std::endl;

    bool mfemHasMPI = false;
  #ifdef MFEM_USE_MPI
    mfemHasMPI = true;
  #endif
    std::cout << "mfem is built " << (mfemHasMPI ? "with" : "without") << " MPI"
              << std::endl;

    if(!axomHasMPI)
    {
      EXPECT_FALSE(mfemHasMPI) << "Axom expects mfem to be built without MPI "
                                  "when it is not built with MPI";
    }
  }

  // Verify that this copy of mfem is configured without Sidre
  {
    bool mfemHasSidre = false;
  #ifdef MFEM_USE_SIDRE
    mfemHasSidre = true;
  #endif

    EXPECT_FALSE(mfemHasSidre) << "Axom expects mfem to be built without Sidre";
  }
}
#endif  // AXOM_USE_MFEM
