#include "lumberjack/Logger.hpp"
#include "lumberjack/RootCommunicator.hpp"
#include "lumberjack/MessageInfo.hpp"

#include "common/CommonTypes.hpp"

#include <mpi.h>
#include <iostream>

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);

	int commRank = -1;
	int commSize = -1;

	MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
	MPI_Comm_rank(MPI_COMM_WORLD, &commSize);

	asctoolkit::lumberjack::RootCommunicator communicator;
	asctoolkit::lumberjack::Logger logger;
	logger.initialize(MPI_COMM_WORLD, &communicator);

	logger.queueMessage("This message is not important");
	logger.pushMessagesOnce();
	MPI_Barrier(MPI_COMM_WORLD);
	std::vector<asctoolkit::lumberjack::MessageInfo>* messageInfos = logger.getMessages();

	if (commRank == 0){
		std::cout << "Rank 0: Printing Messages Recieved!!" << std::endl;
		for(int i=0; i<(int)(messageInfos->size()); ++i){
			std::cout << "   " << messageInfos->at(i).message() << std::endl;
		}
	}
	else {
		if (messageInfos == ATK_NULLPTR){
			std::cout << "Rank " << commRank << ": Success! No messages left in queue!" << std::endl;
		}
		else {
			std::cout << "Rank " << commRank << ": Failure! Messages left in queue!" << std::endl;
		}
	}

	delete messageInfos;

    logger.finalize();
	MPI_Finalize();

	return 0;
}
