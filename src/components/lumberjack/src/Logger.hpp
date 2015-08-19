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
 * \file Logger.hpp
 * \author Chris White (white238@llnl.gov)
 *
 * \brief This file contains the class definition of the Logger. This class
 * is the main class users will interact with.
 *******************************************************************************
 */

#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <string>

#include "mpi.h"
#include "lumberjack/Combiner.hpp"
#include "lumberjack/Communicator.hpp"
#include "lumberjack/MessageEqualityCombiner.hpp"
#include "lumberjack/MessageInfo.hpp"

namespace asctoolkit {
namespace lumberjack {

/*!
 *******************************************************************************
 * \class Logger
 *
 * \brief Class that all user interactions with lumberjack are through.
 *
 *  Concrete instances need to inherit from this class and implement these functions.
 *  You will need to add your Communicator using Logger::Initialize
 *
 * \see BinaryTreeCommunicator RootCommunicator MessageInfo Combiner
 *******************************************************************************
 */
class Logger {
    public:
        /*!
         *****************************************************************************
         * \brief Called to initialize the Logger.
         *
         * This performs any setup work the Logger needs before doing any work.
         * It is required that this is called before using the Logger.
         *
         * \param [in] communicator The lumberjack communicator that will send/recieve messages
         * \param [in] ranksLimit Limit on how many ranks are individually tracker per MessageInfo.
         *****************************************************************************
         */
        void initialize(Communicator* communicator, int ranksLimit);

        /*!
         *****************************************************************************
         * \brief Called to finalize the Logger.
         *
         * This performs any cleanup work the Logger needs to do before going away.
         * It is required that this is the last function called by the Logger.
         *****************************************************************************
         */
        void finalize();

        /*!
         *****************************************************************************
         * \brief Adds a Combiner to the Logger.
         *
         * The Logger can have multiple Combiner's.  This is helpful when you have different
         * combining criteria for different messages. If you try to add the same combiner
         * more than once, the second combiner will not be added.  This is determined solely
         * on Combiner::id. Combiner's that are added first have precedence.  The Logger
         * will call delete on the Combiner when finalize is called.
         *
         * \param [in] combiner The Combiner that will be added.
         *****************************************************************************
         */
        void addCombiner(Combiner* combiner);

        /*!
         *****************************************************************************
         * \brief Removes a Combiner from the Logger.
         *
         * This removes and calls delete on a Combiner held by the Logger. If no 
         * Combiner::id matches the given identifier than nothing is removed. 
         *
         * \param [in] combinerIdentifier The Combiner identifier that will be removed.
         *****************************************************************************
         */
        void removeCombiner(const std::string& combinerIdentifier);

        /*!
         *****************************************************************************
         * \brief Clears all Combiners from the Logger.
         *
         * This removes and calls delete on all Combiners held by the Logger. 
         *****************************************************************************
         */
        void clearCombiners();

        /*!
         *****************************************************************************
         * \brief Fills the given vector with all currently held MessageInfos.
         *
         * This swaps the given vector with this ranks's internal vector that holds it's MessageInfos.
         * The given vector should be empty.  If the rank is not supposed to output Messages,
         * then the vector is not swapped with the internal vector.
         *
         * \param [in,out] filledVector An empty vector that will be filled with MessageInfos.
         *****************************************************************************
         */
        void getMessageInfos(std::vector<MessageInfo*>& filledVector);

        /*!
         *****************************************************************************
         * \brief Sets the rank limit.
         *
         * This is the limit on how many ranks generated a given message are individually tracked
         * per MessageInfo.  After the limit has been reached, only the MessageInfo::rankCount is 
         * incremented.
         *
         * \param [in] ranksLimit Limit on how many ranks are individually tracked per MessageInfo.
         *****************************************************************************
         */
        void ranksLimit(int ranksLimit);

        /*!
         *****************************************************************************
         * \brief Returns the rank limit.
         *
         * This is the limit on how many ranks generated a given message are individually tracked
         * per MessageInfo.  After the limit has been reached, only the MessageInfo::rankCount is 
         * incremented.
         *****************************************************************************
         */
        int ranksLimit();

        /*!
         *****************************************************************************
         * \brief Queues a message to be sent and combined by lumberjack
         *
         * This creates a MessageInfo and queues it to be sent through the Communicator
         * to the root node.  This message may be combined with others depending on the
         * given criteria by the already defined Combiners. Depending on the behavior
         * of the Communicator, the message will not be outputted immediately.
         *
         * \param [in] message Message to be queued.
         *****************************************************************************
         */
        void queueMessage(const std::string& message);

        /*!
         *****************************************************************************
         * \brief Queues a message to be sent and combined by lumberjack
         *
         * This creates a MessageInfo and queues it to be sent through the Communicator
         * to the root node.  This message may be combined with others depending on the
         * given criteria by the already defined Combiners. Depending on the behavior
         * of the Communicator, the message will not be outputted immediately.
         *
         * \param [in] message Message string
         * \param [in] fileName File name of message
         * \param [in] lineNumber Line number of message
         *****************************************************************************
         */
        void queueMessage(const std::string& message, const std::string& fileName, const int lineNumber);

        /*!
         *****************************************************************************
         * \brief This pushes all messages once up the Communicator's tree structure.
         *
         * All of the children push their MessageInfo's to their parent node. This is
         * is helpful if you want to spread the work load of lumberjack over a large
         * set of work. After MessageInfos are pushed they are combined.
         *****************************************************************************
         */
        void pushMessageInfosOnce();

        /*!
         *****************************************************************************
         * \brief This pushes all messages fully up the Communicator's tree structure.
         *
         * All messages are continually pushed until all messages are pushed to the root node.
         * After MessageInfos are pushed once they are combined before being pushed again.
         *****************************************************************************
         */
        void pushMessageInfosFully();
    private:
        /*!
         *****************************************************************************
         * \brief This combines all currently held MessageInfos.
         *
         * All MessageInfos are combined by the currently held Combiners.
         *****************************************************************************
         */
        void combineMessageInfos();

        Communicator* m_communicator;
        int m_ranksLimit;
        MessageEqualityCombiner* m_messageEqualityCombiner;
        std::vector<Combiner*> m_combiners;
        std::vector<MessageInfo*> m_messages;
};

} // end namespace lumberjack
} // end namespace asctoolkit

#endif
