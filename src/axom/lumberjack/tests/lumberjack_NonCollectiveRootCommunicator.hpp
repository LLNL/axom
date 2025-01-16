// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/lumberjack/NonCollectiveRootCommunicator.hpp"
#include "gtest/gtest.h"
#include "mpi.h"

TEST(lumberjack_NonCollectiveRootCommunicator, noncollective_communication)
{
  MPI_Barrier(MPI_COMM_WORLD);

  int commSize = -1;
  MPI_Comm_size(MPI_COMM_WORLD, &commSize);

  const int ranksLimit = 5;
  axom::lumberjack::NonCollectiveRootCommunicator c;

  c.initialize(MPI_COMM_WORLD, ranksLimit);

  std::string message = std::to_string(c.rank());

  std::vector<const char*> receivedPackedMessages;

  // send message only from even ranks that are non-zero
  if((c.rank() % 2) == 0 && c.rank() != 0)
  {
    c.push(message.c_str(), receivedPackedMessages);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  // receive messages from rank 0 after barrier
  if(c.rank() == 0)
  {
    c.push(message.c_str(), receivedPackedMessages);
  }

  if(c.rank() != 0)
  {
    EXPECT_EQ((int)receivedPackedMessages.size(), 0);
  }
  else
  {
    const int numMessagesToReceive =
      ((commSize % 2) == 0) ? ((commSize / 2) - 1) : (commSize / 2);
    EXPECT_EQ((int)receivedPackedMessages.size(), numMessagesToReceive);
    for(int i = 1; i <= numMessagesToReceive; ++i)
    {
      std::string currMessage = std::to_string(i * 2);
      bool found = false;
      for(auto& rm : receivedPackedMessages)
      {
        if(strcmp(rm, currMessage.c_str()) == 0)
        {
          found = true;
        }
      }
      if(!found)
      {
        std::cout << "Error: Message not received: " << currMessage << std::endl;
      }
      EXPECT_EQ(found, true);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  c.finalize();

  MPI_Barrier(MPI_COMM_WORLD);

  // cleanup allocated memory from received messages
  for(auto& rm : receivedPackedMessages)
  {
    delete rm;
  }

  receivedPackedMessages.clear();
}

TEST(lumberjack_NonCollectiveRootCommunicator, multiple_communicators)
{
  MPI_Barrier(MPI_COMM_WORLD);

  int commSize = -1;
  MPI_Comm_size(MPI_COMM_WORLD, &commSize);

  const int ranksLimit = 5;
  axom::lumberjack::NonCollectiveRootCommunicator c1;
  axom::lumberjack::NonCollectiveRootCommunicator c2;

  c1.initialize(MPI_COMM_WORLD, ranksLimit);
  c2.initialize(MPI_COMM_WORLD, ranksLimit);

  // each communicator has unique MPI tag
  EXPECT_NE(c1.mpiTag(), c2.mpiTag());

  MPI_Barrier(MPI_COMM_WORLD);

  c1.finalize();
  c2.finalize();

  MPI_Barrier(MPI_COMM_WORLD);
}
