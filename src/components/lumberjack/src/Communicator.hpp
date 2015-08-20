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
 * \file Communicator.hpp
 * \author Chris White (white238@llnl.gov)
 *
 * \brief This file contains the pure virtual base class definition of the 
 * Communicator.
 *******************************************************************************
 */

#ifndef COMMUNICATOR_HPP
#define COMMUNICATOR_HPP

#include <vector>

#include "lumberjack/Message.hpp"

namespace asctoolkit {
namespace lumberjack {

/*!
 *******************************************************************************
 * \class Communicator
 *
 * \brief Abstract base class defining the interface of all Communicator classes.
 *
 *  Concrete instances need to inherit from this class and implement these functions.
 *  You will need to add your Communicator using Logger::initialize
 *
 * \see BinaryTreeCommunicator RootCommunicator Logger
 *******************************************************************************
 */
class Communicator {
    public:
        /*!
         *****************************************************************************
         * \brief Called to initialize the Communicator.
         *
         * This performs any setup work the Communicator needs before doing any work.
         * It is required that this is called before using the Communicator.
         *
         * \param [in] comm The MPI communicator
         * \param [in] ranksLimit Limit on how many ranks are individually tracked per Message.
         *****************************************************************************
         */
        virtual void initialize(MPI_Comm comm, int ranksLimit) = 0;

        /*!
         *****************************************************************************
         * \brief Called to finalize the Communicator.
         *
         * This performs any cleanup work the Communicator needs to do before going away.
         * It is required that this is the last function called by the Communicator.
         *****************************************************************************
         */
        virtual void finalize() = 0;

        /*!
         *****************************************************************************
         * \brief This pushes all messages once up the Communicator class's tree structure.
         *
         * All of the children push their Message classes to their parent node. This is
         * is helpful if you want to spread the work load of Lumberjack over a large
         * set of work.
         *
         * \param [in,out] messages All of this rank's Message classes.
         *****************************************************************************
         */
        virtual void pushMessagesOnce(std::vector<Message*>& messages) = 0;

        /*!
         *****************************************************************************
         * \brief This pushes all messages fully up the Communicator class's tree structure.
         *
         * All Message classes are continually pushed until all Message classes
         * are pushed to the root node.
         *
         * \param [in,out] messages All of this rank's Message classes.
         *****************************************************************************
         */
        virtual void pushMessagesFully(std::vector<Message*>& messages) = 0;

        /*!
         *****************************************************************************
         * \brief Function used by the Logger to indicate whether this node should be
         * outputting messages. The Communicator class's tree structure dictates this.
         *****************************************************************************
         */
        virtual bool shouldMessagesBeOutputted() = 0;

        /*!
         *****************************************************************************
         * \brief Returns the MPI rank of this node
         *****************************************************************************
         */
        virtual int rank() = 0;
};

} // end namespace lumberjack
} // end namespace asctoolkit

#endif
