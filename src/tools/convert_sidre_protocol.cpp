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
 * \file
 * \brief This file contains a utility to convert a Sidre datastore
 *        from the sidre_hdf5 protocol to another supported protocol.
 *
 * Users must supply a path to a sidre_hdf5 rootfile and base name for
 * the output datastores.  Optional command line arguments include
 * a --protocol option (the default is 'json')
 * and  a --strip option to truncate the array data to at most N elements.
 * The strip option also prepends each array with its original size and a filler
 * entry of 0 for integer arrays or nan for floating point arrays.
 * E.g. if the array had 6 entries [1.01. 2.02, 3.03, 4.04, 5.05, 6.06]
 * and the user passed in --strip 3, the array would be converted to
 * [6, nan, 1.01, 2.02, 3.03].
 *
 * \note The strip option is intended as a temporary solution to truncating
 * a dataset to allow easier debugging.  In the future, we intend to separate
 * the conversion and truncation/display functionality into separate utilities.
 *
 * Usage:
 *    ./convert_sidre_protocol --input path_to_datastore_root_file \
 *                             --output path_to_output_datastore \
 *                             [--protocol a_supported_protocol] \
 *                             [--strip N]
 */
#include "mpi.h"

#include "axom/config.hpp"
#include "axom/Types.hpp"
#include "fmt/format.h"
#include "slic/slic.hpp"
#include "slic/LogStream.hpp"

#ifdef AXOM_USE_LUMBERJACK
  #include "slic/LumberjackStream.hpp"
#else
  #include "slic/GenericOutputStream.hpp"
#endif

#include "sidre/SidreTypes.hpp"
#include "sidre/DataStore.hpp"
#include "sidre/Group.hpp"
#include "sidre/Buffer.hpp"
#include "sidre/View.hpp"

#include "sidre/IOManager.hpp"

#include "slam/SizePolicies.hpp"
#include "slam/OffsetPolicies.hpp"
#include "slam/StridePolicies.hpp"
#include "slam/OrderedSet.hpp"

#include <limits>       // for numeric_limits<int>
#include <cstdlib>      // for atoi


using axom::sidre::DataStore;
using axom::sidre::Group;
using axom::sidre::Buffer;
using axom::sidre::View;
using axom::sidre::IOManager;


typedef axom::sidre::IndexType IndexType;
typedef axom::slam::policies::RuntimeSize<IndexType>   SzPol;
typedef axom::slam::policies::ZeroOffset<IndexType> OffPol;
typedef axom::slam::policies::RuntimeStride<IndexType> StrPol;
typedef axom::slam::OrderedSet<SzPol, OffPol, StrPol> ViewSet;

void setupLogging();
void teardownLogging();

/** Simple structure to hold the parsed command line arguments */
struct CommandLineArguments
{
  static const int NUM_SIDRE_PROTOCOLS = 7;
  static const std::string s_validProtocols[NUM_SIDRE_PROTOCOLS];

  std::string m_inputName;
  std::string m_outputName;
  std::string m_protocol;
  int m_numStripElts;

  CommandLineArguments()
    : m_inputName(""),
    m_outputName(""),
    m_protocol(""),
    m_numStripElts(-1)
  {}

  bool hasInputName() const { return !m_inputName.empty(); }
  bool hasOutputName() const { return !m_outputName.empty(); }
  bool hasOutputProtocol() const { return !m_protocol.empty(); }
  bool shouldStripData() const { return m_numStripElts >= 0; }

  /**  Returns the maximum allowed elements in a view of the output datastore */
  int maxEltsPerView() const
  {
    return shouldStripData()
           ? m_numStripElts
           : std::numeric_limits<int>::max();
  }

  /** Checks whether the input string is a valid sidre protocol */
  static bool isValidProtocol(const std::string& protocol)
  {
    return std::find(s_validProtocols, s_validProtocols + NUM_SIDRE_PROTOCOLS,
                     protocol);
  }

  /** Logs usage information for the utility */
  static void usage()
  {
    fmt::MemoryWriter out;
    out << "Usage ./convert_sidre_protocol <options>";
    out.write("\n\t{:<30}{}", "--help", "Output this message and quit");
    out.write("\n\t{:<30}{}", "--input <file>",
              "(required) Filename of input datastore");
    out.write("\n\t{:<30}{}", "--output <file>",
              "(required) Filename of output datastore");
    out.write("\n\t{:<30}{}", "--strip <N>", "Indicates if data in output file should be "
                                             "stripped (to first N entries) (default: off)");
    out.write("\n\t{:<30}{}", "--protocol <str>",
              "Desired protocol for output datastore (default: json)");

    out.write("\n\n\t{: <40}","Available protocols:");
    for(int i=0 ; i< NUM_SIDRE_PROTOCOLS ; ++i)
    {
      out.write("\n\t  {: <50}", s_validProtocols[i]);
    }

    SLIC_INFO( out.str() );
  }

};

const std::string CommandLineArguments::s_validProtocols[] = {
  "json",
  "sidre_hdf5",
  "sidre_conduit_json",
  "sidre_json",
  "conduit_hdf5",
  "conduit_bin",
  "conduit_json",
};


/** Terminates execution */
void quitProgram(int exitCode = 0)
{
  teardownLogging();
  MPI_Finalize();
  exit(exitCode);
}


/**
 * \brief Utility to parse the command line options
 * \return An instance of the CommandLineArguments struct.
 */
CommandLineArguments parseArguments(int argc, char** argv, int myRank)
{
  CommandLineArguments clargs;

  for(int i=1 ; i< argc ; ++i)
  {
    std::string arg(argv[i]);
    if(arg == "--input")
    {
      clargs.m_inputName = std::string(argv[++i]);
    }
    else if(arg == "--output")
    {
      clargs.m_outputName = std::string(argv[++i]);
    }
    else if(arg == "--protocol")
    {
      clargs.m_protocol = std::string(argv[++i]);
    }
    else if(arg == "--strip")
    {
      clargs.m_numStripElts = std::atoi(argv[++i]);
    }
    else if(arg == "--help" || arg == "-h" )
    {
      if(myRank == 0)
      {
        clargs.usage();
      }
      quitProgram();
    }
  }

  // Input file name is required
  bool isValid = true;;
  if(!clargs.hasInputName())
  {
    SLIC_WARNING("Must supply an input datastore root file.");
    isValid = false;
  }

  if(!clargs.hasOutputName())
  {
    SLIC_WARNING("Must supply a filename for the output datastore.");
    isValid = false;
  }


  // Check that protocol is valid or supply one
  if(!clargs.hasOutputProtocol())
  {
    clargs.m_protocol = CommandLineArguments::s_validProtocols[0];
  }
  else
  {
    if( !clargs.isValidProtocol( clargs.m_protocol ) )
    {
      SLIC_WARNING( clargs.m_protocol << " is not a valid sidre protocol.");
      isValid = false;
    }

  }

  if(!isValid)
  {
    if(myRank == 0)
    {
      clargs.usage();
    }
    quitProgram(1);
  }

  return clargs;
}


/**
 * \brief Helper function to allocate storage for the external data of the input
 * datastore
 *
 * Iterates recursively through the views and groups of the provided group to
 * find the external data views and allocates the required storage within the
 * extPtrs vector
 *
 * \param grp  The group to traverse
 * \param extPtrs [out] A vector to hold pointers to the allocated data
 *
 * \note We also set the data in each allocated array to zeros
 */
void allocateExternalData(Group* grp, std::vector<void*>& extPtrs)
{
  using namespace axom;

  // for each view
  for(sidre::IndexType idx =  grp->getFirstValidViewIndex() ;
      sidre::indexIsValid(idx) ;
      idx = grp->getNextValidViewIndex(idx) )
  {
    View* view = grp->getView(idx);
    if(view->isExternal())
    {
      SLIC_INFO("External view " << view->getPathName()
                                 << " has " << view->getNumElements() << " elements "
                                 << "(" << view->getTotalBytes() << " bytes)."
                );

      const int idx = extPtrs.size();
      const int sz = view->getTotalBytes();
      extPtrs.push_back(new char[sz]);
      std::memset(extPtrs[idx], 0, sz);
      view->setExternalDataPtr( extPtrs[ idx ]);
    }
  }

  // for each group
  for(sidre::IndexType idx =  grp->getFirstValidGroupIndex() ;
      sidre::indexIsValid(idx) ;
      idx = grp->getNextValidGroupIndex(idx) )
  {
    allocateExternalData(grp->getGroup(idx), extPtrs);
  }
}

/**
 * Shifts the data to the right by two elements,
 * The new first value will be the size of the original array
 * The next values will be 0 for integer data and Nan for float data
 * This is followed by the initial values in the original dataset
 *
 * \param view The array view on which we are operating
 * \param origSize The size of the original array
 */
template<typename sidre_type>
void modifyFinalValuesImpl(View* view, int origSize)
{
  SLIC_DEBUG("Looking at view " << view->getPathName());

  sidre_type* arr = view->getData();

  // Uses a Slam set to help manage the indirection to the view data
  // Note: offset is zero since getData() already accounts for the offset
  ViewSet idxSet = ViewSet::SetBuilder()
                   .size(view->getNumElements())
                   .stride(view->getStride());

  #ifdef AXOM_DEBUG
  fmt::MemoryWriter out_fwd;
  for(int i=0 ; i < idxSet.size() ; ++i)
  {
    out_fwd.write("\n\ti: {}; set[i]: {}; arr [ set[i] ] = {}",
                  i, idxSet[i], arr[ idxSet[i] ] );
  }
  SLIC_DEBUG( out_fwd.str() );
  #endif

  // Shift the data over by two
  const int local_offset = 2;
  for(int i=idxSet.size()-1 ; i >= local_offset ; --i)
  {
    arr[ idxSet[i]] = arr[ idxSet[i-local_offset]];
  }

  // Set the first two elements
  arr[ idxSet[0] ] = static_cast<sidre_type>(origSize);
  arr[ idxSet[1] ] = std::numeric_limits<sidre_type>::quiet_NaN();

  #ifdef AXOM_DEBUG
  fmt::MemoryWriter out_rev;
  for(int i=0 ; i < idxSet.size() ; ++i)
  {
    out_rev.write("\n\ti: {}; set[i]: {}; arr [ set[i] ] = {}",
                  i, idxSet[i], arr[ idxSet[i] ] );
  }
  SLIC_DEBUG( out_rev.str() );
  #endif

}


void modifyFinalValues(View* view, int origSize)
{
  SLIC_DEBUG("Truncating view " << view->getPathName());

  using namespace axom;

  switch(view->getTypeID())
  {
  case sidre::INT8_ID:
    modifyFinalValuesImpl<common::int8>(view, origSize);
    break;
  case sidre::INT16_ID:
    modifyFinalValuesImpl<common::int16>(view, origSize);
    break;
  case sidre::INT32_ID:
    modifyFinalValuesImpl<common::int32>(view, origSize);
    break;
  case sidre::INT64_ID:
    modifyFinalValuesImpl<common::int64>(view, origSize);
    break;
  case sidre::UINT8_ID:
    modifyFinalValuesImpl<common::uint8>(view, origSize);
    break;
  case sidre::UINT16_ID:
    modifyFinalValuesImpl<common::uint16>(view, origSize);
    break;
  case sidre::UINT32_ID:
    modifyFinalValuesImpl<common::uint32>(view, origSize);
    break;
  case sidre::UINT64_ID:
    modifyFinalValuesImpl<common::uint64>(view, origSize);
    break;
  case sidre::FLOAT32_ID:
    modifyFinalValuesImpl<common::float32>(view, origSize);
    break;
  case sidre::FLOAT64_ID:
    modifyFinalValuesImpl<common::float64>(view, origSize);
    break;
  default:
    break;
  }
}

/**
 * \brief Recursively traverse views and groups in grp and truncate views to
 * have at most maxSize+2 elements.
 *
 * Within the truncated arrays, the first element will be the size of the
 * original array and the second will
 * be 0 for integers and nan for floating points.
 * This will be followed by (at most) the first maxSize elements of the original
 * array
 */
void truncateBulkData(Group* grp, int maxSize)
{
  using namespace axom;

  // Add two to maxSize
  for(sidre::IndexType idx =  grp->getFirstValidViewIndex() ;
      sidre::indexIsValid(idx) ;
      idx = grp->getNextValidViewIndex(idx) )
  {
    View* view = grp->getView(idx);
    bool isArray = view->hasBuffer() || view->isExternal();

    if(isArray)
    {
      const int numOrigElts = view->getNumElements();
      const int newSize = std::min(maxSize+2, numOrigElts);

      if(view->hasBuffer() && numOrigElts > newSize)
      {
        const int viewStride = view->getStride();
        const int viewOffset = view->getOffset();
        view->apply(newSize,viewOffset, viewStride);
      }
      // external
      else if(view->isExternal() && numOrigElts > newSize)
      {
        void* dataPtr = view->getVoidPtr();
        view->setExternalDataPtr(view->getTypeID(),newSize, dataPtr);
      }

      modifyFinalValues(view, numOrigElts);
    }
  }

  // for each group
  for(sidre::IndexType idx =  grp->getFirstValidGroupIndex() ;
      sidre::indexIsValid(idx) ;
      idx = grp->getNextValidGroupIndex(idx) )
  {
    truncateBulkData(grp->getGroup(idx), maxSize);
  }
}

/** Sets up the logging using lumberjack */
void setupLogging()
{
  using namespace axom;

  slic::initialize();

  slic::setLoggingMsgLevel(slic::message::Info);

#ifdef AXOM_USE_LUMBERJACK
  std::string rankStr = "[<RANK>]";
#else
  std::string rankStr = "";
#endif

  // Formatting for warning, errors and fatal message
  fmt::MemoryWriter wefFmt;
  wefFmt << "\n***********************************\n"
         << rankStr << "[<LEVEL> in line <LINE> of file <FILE>]\n"
         <<"MESSAGE=<MESSAGE>\n"
         << "***********************************\n";
  std::string wefFormatStr = wefFmt.str();

  // Simple formatting for debug and info messages
  std::string diFormatStr = rankStr + "[<LEVEL>]: <MESSAGE>\n";

  slic::LogStream* wefStream;
  slic::LogStream* diStream;

#ifdef AXOM_USE_LUMBERJACK
  const int ranksLimit = 16;
  wefStream = new slic::LumberjackStream( &std::cout, MPI_COMM_WORLD,
                                          ranksLimit, wefFormatStr );
  diStream =  new slic::LumberjackStream( &std::cout, MPI_COMM_WORLD,
                                          ranksLimit, diFormatStr );
#else
  wefStream = new slic::GenericOutputStream( &std::cout, wefFormatStr );
  diStream = new slic::GenericOutputStream( &std::cout, diFormatStr );
#endif

  slic::addStreamToMsgLevel(wefStream, slic::message::Error);
  slic::addStreamToMsgLevel(wefStream, slic::message::Warning);
  slic::addStreamToMsgLevel(diStream,  slic::message::Info);
  slic::addStreamToMsgLevel(diStream,  slic::message::Debug);
}

/** Finalizes logging and flushes streams */
void teardownLogging()
{
  axom::slic::finalize();
}


int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);


  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  setupLogging();

  // parse the command arguments
  CommandLineArguments args = parseArguments(argc, argv, my_rank);

  // Load the original datastore
  SLIC_INFO("Loading datastore from " << args.m_inputName);
  DataStore ds;
  IOManager manager(MPI_COMM_WORLD);
  manager.read(ds.getRoot(), args.m_inputName);
  int num_files = manager.getNumFilesFromRoot(args.m_inputName);

  // Restore any external data pointers
  SLIC_INFO("Loading external data from datastore");
  std::vector<void*> externalDataPointers;
  allocateExternalData(ds.getRoot(), externalDataPointers);
  manager.loadExternalData(ds.getRoot(), args.m_inputName);

  axom::slic::flushStreams();

  // Internal processing
  if(args.shouldStripData())
  {
    const int numElts = args.maxEltsPerView();
    SLIC_INFO("Truncating views to at most " << numElts << " elements.");

    truncateBulkData(ds.getRoot(), numElts);

    // Add a string view to the datastore to indicate that we modified the data
    fmt::MemoryWriter fout;
    fout << "This datastore was created by the convert_sidre_protocol "
         << "utility with option '--strip " << numElts
         << "'. To simplify debugging, the bulk data in this datastore "
         << "has been truncated to have at most " << numElts
         << " actual values. The first two values are size of the original "
         << " array followed by a zero/Nan. These are followed by (at most) "
         << "the first " << numElts << " values of the array.";

    ds.getRoot()->createViewString("Note", fout.str());
  }

  // Write out datastore to the output file in the specified protocol
  SLIC_INFO("Writing out datastore in "
            << args.m_protocol << " protocol to file(s) with base name "
            << args.m_outputName);
  manager.write(ds.getRoot(), num_files, args.m_outputName, args.m_protocol);


  // Free up memory associated with external data
  for(std::vector<void*>::iterator it = externalDataPointers.begin() ;
      it != externalDataPointers.end() ;
      ++it)
  {
    delete [] static_cast<char*>(*it);
    *it = AXOM_NULLPTR;
  }

  teardownLogging();
  MPI_Finalize();
  return 0;
}
