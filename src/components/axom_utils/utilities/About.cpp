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

#include "axom_utils/utilities/About.hpp"
#include "axom/config.hpp"

#include <algorithm>  // copy
#include <iostream>   // for cout
#include <sstream>    // for stringstream
#include <string>     // for C++ string
#include <vector>     // for STL vector
#include <iterator>   // for ostream_iterator


namespace axom
{

//-----------------------------------------------------------------------------
void about()
{
  about(std::cout);
}

//-----------------------------------------------------------------------------
// The guts of this were taken from the utils_config unit test.
//-----------------------------------------------------------------------------
void about(std::ostream &oss)
{
  oss << "HERE'S SOME INTERESTING STUFF:" << std::endl << std::endl;

  oss << "Axom (" << AXOM_VERSION_FULL << ")"  << std::endl << std::endl;

  oss << "AXOM_VERSION_FULL: "  << AXOM_VERSION_FULL  << std::endl;
  oss << "AXOM_VERSION_MAJOR: " << AXOM_VERSION_MAJOR << std::endl;
  oss << "AXOM_VERSION_MINOR: " << AXOM_VERSION_MINOR << std::endl;
  oss << "AXOM_VERSION_PATCH: " << AXOM_VERSION_PATCH << std::endl;

#ifdef AXOM_VERSION_EXTRA
  oss << "AXOM_VERSION_EXTRA: " << AXOM_VERSION_EXTRA << std::endl;
#endif

  oss << "Compiler Settings: "  << std::endl
      << "   C++11 support: "
#ifdef AXOM_USE_CXX11
    << "ENABLED"
#else
    << "DISABLED"
#endif
    << std::endl
    << "   OpenMP support: "
#ifdef AXOM_USE_OPENMP
    << "ENABLED"
#else
    << "DISABLED"
#endif
    << std::endl;

  oss << "Available components: " << std::endl;

  std::vector<std::string> comps;

#ifdef AXOM_USE_AXOM_UTILS
  comps.push_back("axom_utils");
#endif

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
  std::copy( comps.begin(), comps.end(),
             std::ostream_iterator<std::string>(sstr, "; "));
  oss << " { " << sstr.str() << "}" << std::endl;

  oss << "Active Dependencies: " << std::endl;

  std::vector<std::string> libs;

#ifdef AXOM_USE_CONDUIT
  libs.push_back("conduit");
#endif

#ifdef AXOM_USE_FMT
  libs.push_back("fmt");
#endif

#ifdef AXOM_USE_HDF5
  libs.push_back("hdf5");
#endif

#ifdef AXOM_USE_MFEM
  libs.push_back("mfem");
#endif

#ifdef AXOM_USE_MPI
  libs.push_back("mpi");
#endif

#ifdef AXOM_USE_SPARSEHASH
  libs.push_back("sparsehash");
#endif

  // reset sstr
  sstr.str("");
  std::copy( libs.begin(), libs.end(),
             std::ostream_iterator<std::string>(sstr, "; "));
  oss << " { " << sstr.str() << "}" << std::endl;
}


} // end namespace axom
