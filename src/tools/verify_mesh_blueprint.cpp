// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)


/**
 * \file verify_mesh_blueprint.cpp
 * \brief This file contains a utility to verify that a Sidre datastore
 *        contains a valid mesh blueprint at a given path.
 */

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/sidre.hpp"
#include "axom/slam.hpp"

#include "conduit_blueprint.hpp"

#include "fmt/fmt.hpp"
#include "CLI11/CLI11.hpp"

#include <sstream>
#include <fstream>
#include <utility>
#include <string>

namespace sidre = axom::sidre;
namespace slam = axom::slam;
namespace slic = axom::slic;

void setupLogging();
void teardownLogging();

/** Simple structure to hold the parsed command line arguments */
struct CommandLineArguments
{
  static const std::set<std::string> s_validProtocols;

  std::string m_inputName;
  std::string m_outputPrefix;
  std::string m_blueprintPath;
  bool m_shouldExplore;

  CommandLineArguments()
    : m_inputName("")
    , m_outputPrefix("")
    , m_blueprintPath("")
    , m_shouldExplore(false)
  {}

  void parse(int argc, char** argv, CLI::App& app);

  bool shouldExplore() const { return m_shouldExplore; }
};


/** Terminates execution */
void quitProgram(int exitCode = 0)
{
  teardownLogging();
  MPI_Finalize();
  exit(exitCode);
}


/** Parse the command line arguments */
void CommandLineArguments::parse(int argc, char** argv, CLI::App& app)
{
  app.add_option("-i,--input", m_inputName,
                 "Rootfile for sidre datastore")
  ->required()
  ->check(CLI::ExistingFile);

  app.add_option("-o,--output", m_outputPrefix,
                 "Prefix for dumping blueprint verification."
                 "Will skip if not present.");

  app.add_option("-p,--path", m_blueprintPath,
                 "Path to root of blueprint");
  // TODO: check that this is a valid path? E.g. characters and '/'.

  app.add_flag("-e,--explore", m_shouldExplore,
               "List child groups and views from provided path and quit");

  bool verboseOutput = false;
  app.add_flag("-v,--verbose", verboseOutput,
               "Sets output to verbose")
  ->capture_default_str();

  app.get_formatter()->column_width(35);

  // Could throw an exception
  app.parse(argc, argv);

  if(verboseOutput)
  {
    slic::setLoggingMsgLevel(slic::message::Debug);
  }
}

/** Sets up the logging using lumberjack */
void setupLogging()
{
  slic::initialize();
  slic::setLoggingMsgLevel(slic::message::Info);

#ifdef AXOM_USE_LUMBERJACK
  std::string rankStr = "[<RANK>]";
#else
  std::string rankStr = "";
#endif

  // Formatting for warning, errors and fatal message
  std::stringstream wefFmt;
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

  // the following is helpful for debugging
  // slic::debug::checksAreErrors = true;
}

/** Finalizes logging and flushes streams */
void teardownLogging()
{
  slic::finalize();
}

/// Returns true if \a val is true on all ranks
bool allTrue(bool val)
{
  int localVal = val ? 0 : 1;
  int globalVal = 0;

  MPI_Allreduce(&localVal, &globalVal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  return globalVal == 0;
}

/// Returns true if \a val is true on any rank
bool anyTrue(bool val)
{
  int localVal = val ? 1 : 0;
  int globalVal = 0;

  MPI_Allreduce(&localVal, &globalVal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  return globalVal > 0;
}


int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  setupLogging();


  // parse the command line arguments
  CommandLineArguments args;
  CLI::App app {"Sidre protocol converter"};

  try
  {
    args.parse(argc, argv, app);
  }
  catch (const CLI::ParseError &e)
  {
    int retval = -1;
    if(my_rank==0)
    {
      retval = app.exit(e);
    }
    MPI_Bcast(&retval, 1, MPI_INT, 0, MPI_COMM_WORLD);
    quitProgram(retval);
  }


  // Load the datastore
  SLIC_INFO("Loading datastore from " << args.m_inputName);
  sidre::DataStore ds;
  sidre::IOManager manager(MPI_COMM_WORLD);
  manager.read(ds.getRoot(), args.m_inputName);

  int num_files = manager.getNumFilesFromRoot(args.m_inputName);
  SLIC_INFO("Datastore is stored in " << num_files << " files.");

  // Get the requested group
  bool validPath = false;
  auto* root = ds.getRoot();
  auto* grp = args.m_blueprintPath.empty()
              ? root
              : root->getGroup(args.m_blueprintPath);
  if(grp == nullptr)
  {
    SLIC_WARNING("Invalid path: '" << args.m_blueprintPath << "'.");
  }
  else
  {
    validPath = true;
  }

  bool anyRankValid = anyTrue(validPath);
  if (!anyRankValid )
  {
    SLIC_WARNING("Path was not valid on any rank");
    quitProgram(1);
  }


  // If in explore mode, print out names of groups and views rooted at grp
  if(args.shouldExplore())
  {
    SLIC_INFO("Exploring datastore at path '" << args.m_blueprintPath << "'.");

    std::stringstream sstr;

    if(validPath)
    {
      // loop through groups
      {
        sstr << "Groups:" << std::endl;
        for(auto idx =  grp->getFirstValidGroupIndex() ;
            sidre::indexIsValid(idx) ;
            idx = grp->getNextValidGroupIndex(idx) )
        {
          sstr << "\t" << grp->getGroup(idx)->getPathName() << std::endl;
        }
      }

      // loop through views
      {
        sstr << "Views:" << std::endl;
        for(auto idx =  grp->getFirstValidViewIndex() ;
            sidre::indexIsValid(idx) ;
            idx = grp->getNextValidViewIndex(idx) )
        {
          auto* v  = grp->getView(idx);

          sstr << "\t" << v->getPathName() << " -- ";
          v->print(sstr);
          sstr << std::endl;
        }
      }
      SLIC_INFO("\n" << sstr.str());
    }

    quitProgram();
  }

  // otherwise, verify that the path corresponds to a blueprint.
  bool success = false;
  {
    conduit::Node mesh_node;
    grp->createNativeLayout(mesh_node);

    conduit::Node info;
    if(validPath)
    {
      if(conduit::blueprint::verify("mesh", mesh_node, info))
      {
        SLIC_INFO("Valid blueprint at path '" << args.m_blueprintPath << "'.");
        success = true;
      }
      else
      {
        SLIC_INFO("Invalid blueprint at path '" << args.m_blueprintPath <<
        "'.");
        success = false;
      }

      if(!args.m_outputPrefix.empty() )
      {
        int myRank;
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

        std::string fname = fmt::format("{}_{:04}.log", args.m_outputPrefix,
                                        myRank);

        std::fstream fstr(fname, std::ios::out);
        if(!fstr.is_open())
        {
          SLIC_WARNING("Could not open '" << fname << "'");
        }
        else
        {
          SLIC_INFO("CWD: " << axom::utilities::filesystem::getCWD() );
          SLIC_INFO( "Writing blueprint verify info to '" << fname << "'");
          info.to_json_stream(fstr);
        }

      }
    }

    // MPI reduce the success variable;
    success = allTrue(success);
  }

  teardownLogging();
  MPI_Finalize();
  return success ? 0 : 1;
}
