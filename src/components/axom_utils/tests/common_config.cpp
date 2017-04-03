/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


#include "gtest/gtest.h"

#include "axom/config.hpp"

#ifdef AXOM_USE_OPENMP
  #include <omp.h>
#endif

#ifdef AXOM_USE_BOOST
  #include "boost/version.hpp"
#endif

#include <string>
#include <sstream>          // stringstream
#include <iostream>         // cout
#include <vector>
#include <algorithm>        // copy


TEST(gtest_common_config,config_libraries)
{
  // This test checks which libraries are available in the configuration

  std::cout << "Available libraries: " << std::endl;

  std::vector<std::string> libs;

#ifdef AXOM_USE_BOOST
  libs.push_back("boost");
#endif

#ifdef AXOM_USE_CONDUIT
  libs.push_back("conduit");
#endif

#ifdef AXOM_USE_CXX11
  libs.push_back("C++11");
#endif

#ifdef AXOM_USE_FMT
  libs.push_back("fmt");
#endif

#ifdef AXOM_USE_HDF5
  libs.push_back("hdf5");
#endif

#ifdef AXOM_USE_OPENMP
  libs.push_back("openmp");
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

#ifdef AXOM_USE_STD_UNORDERED_MAP
  libs.push_back("std::unordered_map");
#endif

  std::stringstream sstr;
  std::copy( libs.begin(), libs.end(), std::ostream_iterator<std::string>(sstr, "; "));
  std::cout << "\t{ " << sstr.str() << "}" << std::endl;

  EXPECT_TRUE(true);
}


TEST(gtest_common_config,config_components)
{
  // This test checks which toolkit components are available in the configuration

  std::cout << "Available components: " << std::endl;

  std::vector<std::string> comps;

#ifdef AXOM_USE_COMMON
  comps.push_back("common");
#endif

  EXPECT_EQ(1u, comps.size()) << "Common component is always available.";

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

#ifdef AXOM_USE_SPIO
  comps.push_back("spio");
#endif


  std::stringstream sstr;
  std:: copy( comps.begin(), comps.end(), std::ostream_iterator<std::string>(sstr, "; "));
  std::cout << "\t{ " << sstr.str() << "}" << std::endl;

  EXPECT_TRUE(true);
}

TEST(gtest_common_config,config_openmp)
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
      std::cout <<"\tMy thread id is: " << thId
              <<"\tNum threads is: " << thNum
              <<"\tMax threads is: " << thMax
              << std::endl;
    }

#else
    std::cout << "OpenMP not available in this configuration." << std::endl;
#endif


    // Tests a simple reduction over integers
    // Note: Non-openmp configurations are compiled with a flag
    //       to ignore usage of unknown openmp pragmas
    const int N = 100;
    int sum=0;
    #pragma omp parallel for reduction(+:sum)
    for(int i=1; i<=N; ++i)
    {
        sum += i;
    }
    EXPECT_EQ( N * (N+1) / 2, sum) << "Bad reduction on first " << N << " integers";
    std::cout << "Sum of first " << N << " numbers is: " << sum << std::endl;

}

#ifdef AXOM_USE_BOOST
TEST(gtest_common_config,boost_version)
{
    std::cout << "Using boost version "
          << BOOST_VERSION / 100000     << "."  // major version
          << BOOST_VERSION / 100 % 1000 << "."  // minor version
          << BOOST_VERSION % 100                // patch level
          << std::endl;
}
#endif // AXOM_USE_BOOST

