#ifndef ROOTCOMMUNICATOR_HPP
#define ROOTCOMMUNICATOR_HPP

#include <string>

#include "mpi.h"

#include "common/CommonTypes.hpp"

#include "lumberjack/Communicator.hpp"
#include "lumberjack/MessageInfo.hpp"

namespace asctoolkit {
namespace lumberjack {

class RootCommunicator: public Communicator {
    public:
        void initialize(MPI_Comm comm, int ranksLimit);
        void finalize();

        void pushMessagesOnce(std::vector<MessageInfo*>& messages);
        void pushMessagesFully(std::vector<MessageInfo*>& messages);

        bool shouldMessagesBeOutputted();
        int rank();
    private:
        MPI_Comm m_mpiComm;
        int m_mpiCommRank;
        int m_mpiCommSize;
        int m_ranksLimit;
};

}
}

#endif