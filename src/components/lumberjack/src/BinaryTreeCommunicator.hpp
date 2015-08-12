/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

#ifndef BINARYTREECOMMUNICATOR_HPP
#define BINARYTREECOMMUNICATOR_HPP

#include <string>

#include "mpi.h"

namespace asctoolkit {
namespace lumberjack {

class BinaryTreeCommunicator: public Communicator {
    public:
        void setCommunicator(MPI_Comm comm);
        void pushMessagesOnceUpTree();
        void pushMessagesFullyUpTree();
        std::vector<MessageInfos> getMessages();
        void queueMessage(const std::string& message,
                          const std::string& fileName,
                          const int lineNumber);
        void queueMessage(const std::string& message);
    private:
        MPI_Comm m_comm;
        std::vector<int> m_ranks;
};

}
}

#endif