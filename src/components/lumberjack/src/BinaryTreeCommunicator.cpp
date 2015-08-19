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
 * \file BinaryTreeCommunicator.cpp
 * \author Chris White (white238@llnl.gov)
 *
 * \brief This file contains the class implementation of the BinaryTreeCommunicator.
 *******************************************************************************
 */

#include "lumberjack/BinaryTreeCommunicator.hpp"

namespace asctoolkit {
namespace lumberjack {

void BinaryTreeCommunicator::setCommunicator(MPI_Comm comm)
{

}

void BinaryTreeCommunicator::pushMessagesOnce()
{

}

void BinaryTreeCommunicator::pushMessagesFully()
{

}

std::vector<MessageInfos> BinaryTreeCommunicator::getMessages()
{

}

void BinaryTreeCommunicator::queueMessage(const std::string& message,
                                    const std::string& fileName,
                                    const int lineNumber)
{
    MessageInfo mi(message, m_commRank, fileName, lineNumber);
    m_messages.append(mi);
}

void BinaryTreeCommunicator::queueMessage(const std::string& message)
{
    queueMessage(message, "", -1);
}

} // end namespace lumberjack
} // end namespace asctoolkit
