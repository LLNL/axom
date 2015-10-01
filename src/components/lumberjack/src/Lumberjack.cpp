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
 * \file Lumberjack.cpp
 * \author Chris White (white238@llnl.gov)
 *
 * \brief This file contains the class implementation of the Lumberjack. This class
 * is the main class users will interact with.
 *******************************************************************************
 */

#include "lumberjack/Lumberjack.hpp"

#include "common/CommonTypes.hpp"

namespace asctoolkit {
namespace lumberjack {

void Lumberjack::initialize(Communicator* communicator, int ranksLimit)
{
    m_communicator = communicator;
    m_ranksLimit = ranksLimit;
    m_combiners.push_back(new TextEqualityCombiner);
}

void Lumberjack::finalize()
{
    m_communicator = ATK_NULLPTR;
    clearCombiners();
}

void Lumberjack::addCombiner(Combiner* combiner)
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

void Lumberjack::removeCombiner(const std::string& combinerIdentifier)
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

void Lumberjack::clearCombiners()
{
    for (int i=0; i<(int)m_combiners.size(); ++i){
        delete m_combiners[i];
    }
    m_combiners.clear();
}

void Lumberjack::getMessages(std::vector<Message*>& filledVector)
{
    if (m_communicator->shouldMessagesBeOutputted()){
        filledVector.swap(m_messages);
    }
}

void Lumberjack::ranksLimit(int value)
{
    m_ranksLimit = value;
    m_communicator->ranksLimit(value);
}

int Lumberjack::ranksLimit()
{
    return m_ranksLimit;
}

void Lumberjack::queueMessage(const std::string& text)
{
    queueMessage(text, "", -1);
}

void Lumberjack::queueMessage(const std::string& text, const std::string& fileName, const int lineNumber)
{
    Message* mi = new Message(text, m_communicator->rank(), fileName, lineNumber);
    m_messages.push_back(mi);
}

void Lumberjack::pushMessagesOnce()
{
    m_communicator->pushMessagesOnce(m_messages, m_combiners);
}

void Lumberjack::pushMessagesFully()
{
    m_communicator->pushMessagesFully(m_messages, m_combiners);
}

} // end namespace lumberjack
} // end namespace asctoolkit
