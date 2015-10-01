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
 * \file BinaryTreeCommunicator.hpp
 * \author Chris White (white238@llnl.gov)
 *
 * \brief This file contains the class definition of the BinaryTreeCommunicator.
 *******************************************************************************
 */

#ifndef BINARYTREECOMMUNICATOR_HPP
#define BINARYTREECOMMUNICATOR_HPP

#include <string>

#include "mpi.h"

#include "lumberjack/Communicator.hpp"
#include "lumberjack/Message.hpp"

namespace asctoolkit {
namespace lumberjack {

/*!
 *******************************************************************************
 * \class BinaryTreeCommunicator
 *
 * \brief A Communicator that utilizes a binary tree stucture of MPI nodes 
 * to scalably pass Message classes.
 *
 *  You will need to add your Communicator using Lumberjack::initialize.
 *
 * \see Communicator Lumberjack
 *******************************************************************************
 */
class BinaryTreeCommunicator: public Communicator {
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
         * \brief All children push their Message classes pushed to their parent and then combined.
         *
         * This is helpful when you want to spread your Lumberjack work over a set of 
         * your work. Instead of doing all of Lumberjack's work at one point in your program.
         * This does not guarantee your Message will be ready to output.  After this call
         * only Message classes at the root node will be returned from Lumberjack::getMessages.
         *
         * \param [in,out] messages All of this rank's Message classes.
         * \param [in,out] combiners All of currently active Combiner classes.
         *****************************************************************************
         */
        void pushMessagesOnce(std::vector<Message*>& messages, std::vector<Combiner*>& combiners);

        /*!
         *****************************************************************************
         * \brief All Message classes pushed from child to parent until all have passed
         * through Lumberjack. Messages are combined after each iteration.
         *
         * This is helpful when you want to completely clear out Message classes in Lumberjack.
         *
         * \param [in,out] messages All of this rank's Message classes.
         * \param [in,out] combiners All of currently active Combiner classes.
         *****************************************************************************
         */
        void pushMessagesFully(std::vector<Message*>& messages, std::vector<Combiner*>& combiners);

        /*!
         *****************************************************************************
         * \brief Function used by the Lumberjack to indicate whether this node should be
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
        int m_treeHeight;
        int m_parentRank;
        int m_leftChildRank;
        int m_rightChildRank;
        int m_childCount;
};

} // end namespace lumberjack
} // end namespace asctoolkit

#endif
