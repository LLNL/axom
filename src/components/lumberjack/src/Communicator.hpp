/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

#ifndef COMMUNICATOR_HPP
#define COMMUNICATOR_HPP

#include <vector>

#include "lumberjack/MessageInfo.hpp"

namespace asctoolkit {
namespace lumberjack {

class Communicator {
    public:
        virtual void initialize(MPI_Comm comm, int ranksLimit) = 0;
        virtual void finalize() = 0;
        virtual void pushMessagesOnce(std::vector<MessageInfo*>& messages) = 0;
        virtual void pushMessagesFully(std::vector<MessageInfo*>& messages) = 0;
        virtual bool shouldMessagesBeOutputted() = 0;
        virtual int rank() = 0;
};

}
}

#endif