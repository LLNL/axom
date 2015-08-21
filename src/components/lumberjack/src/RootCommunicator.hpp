/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

/*!
 *******************************************************************************
 * \file RootCommunicator.hpp
 * \author Chris White (white238@llnl.gov)
 *
 * \brief This file contains the class definition of the RootCommunicator.
 *******************************************************************************
 */

#ifndef ROOTCOMMUNICATOR_HPP
#define ROOTCOMMUNICATOR_HPP

#include <string>

#include "mpi.h"

#include "lumberjack/Communicator.hpp"
#include "lumberjack/Message.hpp"

namespace asctoolkit {
namespace lumberjack {

/*!
 *******************************************************************************
 * \class RootCommunicator
 *
 * \brief A Communicator that all MPI nodes communicate with the root node. This class
 * does NOT scale and is provided for demonstration purposes only.
 *
 *  You will need to add your Communicator using Logger::initialize.
 *
 * \see BinaryTreeCommunicator Communicator Logger
 *******************************************************************************
 */
class RootCommunicator: public Communicator {
    public:
        /*!
         *****************************************************************************
         * \brief Called to initialize the Communicator.
         *
         * This performs any setup work the Communicator needs before doing any work.
         * It is required that this is called before using the Communicator.
         *
         * \param [in] comm The MPI Communicator
         * \param [in] ranksLimit Limit on how many ranks are individually tracked per Message.
         *****************************************************************************
         */
        void initialize(MPI_Comm comm, int ranksLimit);

        /*!
         *****************************************************************************
         * \brief Called to finalize the Communicator.
         *
         * This performs any cleanup work the Communicator needs to do before going away.
         * It is required that this is the last function called by the Communicator.
         *****************************************************************************
         */
        void finalize();

        /*!
         *****************************************************************************
         * \brief Sets the rank limit.
         *
         * This is the limit on how many ranks generated a given message are individually tracked
         * per Message.  After the limit has been reached, only the Message::rankCount is 
         * incremented.
         *
         * \param [in] value Limit on how many ranks are individually tracked per Message.
         *****************************************************************************
         */
        void ranksLimit(int value);

        /*!
         *****************************************************************************
         * \brief Returns the rank limit.
         *
         * This is the limit on how many ranks generated a given message are individually tracked
         * per Message.  After the limit has been reached, only the Message::rankCount is 
         * incremented.
         *****************************************************************************
         */
        int ranksLimit();

        /*!
         *****************************************************************************
         * \brief This pushes all messages to the root node.
         *
         * All messages are pushed to the root node. This is the same as 
         * RootCommunicator::pushMessagesFully for this Communicator.
         *
         * \param [in,out] messages All of this rank's Message classes.
         *****************************************************************************
         */
        void pushMessagesOnce(std::vector<Message*>& messages);

        /*!
         *****************************************************************************
         * \brief This pushes all messages to the root node.
         *
         * All messages are pushed to the root node. This is the same as 
         * RootCommunicator::pushMessagesOnce for this Communicator.
         *
         * \param [in,out] messages All of this rank's Message classes.
         *****************************************************************************
         */
        void pushMessagesFully(std::vector<Message*>& messages);

        /*!
         *****************************************************************************
         * \brief Function used by the Logger to indicate whether this node should be
         * outputting messages. Only the root node outputs messages.
         *****************************************************************************
         */
        bool shouldMessagesBeOutputted();

        /*!
         *****************************************************************************
         * \brief Returns the MPI rank of this node
         *****************************************************************************
         */
        int rank();
    private:
        MPI_Comm m_mpiComm;
        int m_mpiCommRank;
        int m_mpiCommSize;
        int m_ranksLimit;
};

} // end namespace lumberjack
} // end namespace asctoolkit

#endif
