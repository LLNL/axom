/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

/*!
 *******************************************************************************
 * \file parallel_logging_example.cc
 *
 * \date May 7, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *
 *******************************************************************************
 */

// C/C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// Toolkit includes
#include "common/CommonTypes.hpp"
#include "lumberjack/Lumberjack.hpp"
#include "lumberjack/BinaryTreeCommunicator.hpp" 
#include "slic/slic.hpp"
#include "slic/LumberjackStream.hpp"

// MPI
#include <mpi.h>

using namespace asctoolkit;

#define CYCLELIMIT 10
#define RANKSLIMIT 5

//------------------------------------------------------------------------------
int main( int argc, char** argv )
{
  // Process command line options
  if (argc != 2) {
    std::cout << "Error: Wrong amount of command line arguments given. Usage:" << std::endl << 
                 "   " << argv[0] << " <file to be read>" << std::endl;
    return 1;
  }
  char* fileName = argv[1];

  // Initialize MPI
  MPI_Init( &argc, &argv );
  int rank=-1;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // Initialize Lumberjack
  lumberjack::BinaryTreeCommunicator ljComm;
  ljComm.initialize(MPI_COMM_WORLD, RANKSLIMIT);
  asctoolkit::lumberjack::Lumberjack lj;
  lj.initialize(&ljComm, RANKSLIMIT);

  // Initialize SLIC
  std::string format = std::string( "<MESSAGE>\n") +
                       std::string( "\t<TIMESTAMP>" ) +
                       std::string( "\tLEVEL=<LEVEL>\n") +
                       std::string( "\tFILE=<FILE>\n") +
                       std::string( "\tLINE=<LINE>\n");
  slic::initialize();

  // Open file for outputting messages
  std::ofstream outFile;
  outFile.open("speedTestOutput");

  // Set SLIC logging level and stream
  slic::setLoggingMsgLevel( slic::message::Debug );
  slic::LumberjackStream* ljStream = 
        new slic::LumberjackStream( &outFile, &lj, format );
  slic::addStreamToAllMsgLevels( ljStream );

  // Read lines from file
  std::string currMessage;
  std::vector<std::string> lines;
  std::ifstream file(fileName);
  while(std::getline(file, currMessage)){
    currMessage += '\n';
    lines.push_back(currMessage);
  }
  file.close();

  // Queue messages
  int cycleCount = 0;
  int linesSize = (int)lines.size();
  for (int i = 0; i < linesSize; ++i){
    slic::logMessage( slic::message::Info,
                         lines[i],
                         __FILE__,
                         __LINE__
                         );
    ++cycleCount;
    if (cycleCount > CYCLELIMIT) {
      // Incrementally push messages through the log stream
      slic::pushStreams();
      cycleCount = 0;
    }
  }

  // Fully flush system of messages
  slic::flushStreams();

  // Shutdown SLIC and Lumberjack
  slic::finalize();
  ljComm.finalize();
  lj.finalize();

  // Finalize MPI
  MPI_Finalize();

  return 0;
}
