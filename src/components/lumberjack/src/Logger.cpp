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
 * \file Logger.cpp
 * \author Chris White (white238@llnl.gov)
 *
 * \brief This file contains the class implementation of the Logger. This class
 * is the main class users will interact with.
 *******************************************************************************
 */

#include "lumberjack/Logger.hpp"

#include "common/CommonTypes.hpp"

namespace asctoolkit {
namespace lumberjack {

void Logger::initialize(Communicator* communicator, int ranksLimit)
{
    m_communicator = communicator;
    m_ranksLimit = ranksLimit;
    m_combiners.push_back(new TextEqualityCombiner);
}

void Logger::finalize()
{
    m_communicator = ATK_NULLPTR;
    clearCombiners();
}

void Logger::addCombiner(Combiner* combiner)
{
    bool identifierFound = false;
    for (int i=0; i<(int)m_combiners.size(); ++i){
        if (m_combiners[i]->id() == combiner->id()){
            identifierFound = true;
            break;
        }
    }
    if (!identifierFound){
        m_combiners.push_back(combiner);
    }
}

void Logger::removeCombiner(const std::string& combinerIdentifier)
{
    int combinerToBeRemoved = -1;
    for (int i=0; i<(int)m_combiners.size(); ++i){
        if (m_combiners[i]->id() == combinerIdentifier){
            delete m_combiners[i];
            combinerToBeRemoved = i;
            break;
        }
    }
    if (combinerToBeRemoved != -1){
        m_combiners.erase(m_combiners.begin()+combinerToBeRemoved);
    }
}

void Logger::clearCombiners()
{
    for (int i=0; i<(int)m_combiners.size(); ++i){
        delete m_combiners[i];
    }
    m_combiners.clear();
}

void Logger::getMessages(std::vector<Message*>& filledVector)
{
    if (m_communicator->shouldMessagesBeOutputted()){
        filledVector.swap(m_messages);
    }
}

void Logger::ranksLimit(int value)
{
    m_ranksLimit = value;
    m_communicator->ranksLimit(value);
}

int Logger::ranksLimit()
{
    return m_ranksLimit;
}

void Logger::queueMessage(const std::string& text)
{
    queueMessage(text, "", -1);
}

void Logger::queueMessage(const std::string& text, const std::string& fileName, const int lineNumber)
{
    Message* mi = new Message(text, m_communicator->rank(), fileName, lineNumber);
    m_messages.push_back(mi);
}

void Logger::pushMessagesOnce()
{
    m_communicator->pushMessagesOnce(m_messages);
    combineMessages();
}

void Logger::pushMessagesFully()
{
    m_communicator->pushMessagesFully(m_messages);
    combineMessages();
}

void Logger::combineMessages()
{
    int messagesSize = (int)m_messages.size();
    if (messagesSize < 2){
        return;
    }

    std::vector<Message*> finalMessages;
    std::vector<int> indexesToBeDeleted;
    int combinersSize = (int)m_combiners.size();
    bool combinedMessage = false;
    finalMessages.push_back(m_messages[0]);
    for (int allIndex=1; allIndex<messagesSize; ++allIndex){
        combinedMessage = false;
        for (int finalIndex=0; finalIndex<(int)finalMessages.size(); ++finalIndex){
            for (int combinerIndex=0; combinerIndex<combinersSize; ++combinerIndex){
                if (m_combiners[combinerIndex]->shouldMessagesBeCombined(*finalMessages[finalIndex],
                                                                             *m_messages[allIndex])){
                    m_combiners[combinerIndex]->combine(*finalMessages[finalIndex],
                                                        *m_messages[allIndex], m_ranksLimit);
                    indexesToBeDeleted.push_back(allIndex);
                    combinedMessage = true;
                    break;
                }
            }
            if (combinedMessage){
                break;
            }
        }
        if (!combinedMessage){
            finalMessages.push_back(m_messages[allIndex]);
        }
    }

    for (int i=0; i<(int)indexesToBeDeleted.size(); ++i){
        delete m_messages[indexesToBeDeleted[i]];
    }
    m_messages.swap(finalMessages);
}

} // end namespace lumberjack
} // end namespace asctoolkit
