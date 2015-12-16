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

#include "lumberjack/Utility.hpp"

#include <cstring>

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
    clearMessages();
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

void Lumberjack::clearMessages()
{
    for (int i=0; i<(int)m_messages.size(); ++i) {
        delete m_messages[i];
    }
    m_messages.clear();
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
    const char* packedMessagesToBeSent = "";
    if (!m_communicator->shouldMessagesBeOutputted()) {
        combineMessages();
        packedMessagesToBeSent = packMessages();
        clearMessages();
    }
    std::vector<const char*> receivedPackedMessages;

    m_communicator->push(packedMessagesToBeSent, receivedPackedMessages);

    for (int i=0;i<(int)receivedPackedMessages.size(); ++i){
        unpackMessages(receivedPackedMessages[i]);
        delete receivedPackedMessages[i];
    }
    receivedPackedMessages.clear();

    combineMessages();
}

void Lumberjack::pushMessagesFully()
{
    const char* packedMessagesToBeSent = "";
    std::vector<const char*> receivedPackedMessages;
    int numPushesToFlush = m_communicator->numPushesToFlush();
    for (int i=0; i<numPushesToFlush; ++i){
        if (!m_communicator->shouldMessagesBeOutputted()) {
            combineMessages();
            packedMessagesToBeSent = packMessages();
            clearMessages();
        }

        m_communicator->push(packedMessagesToBeSent, receivedPackedMessages);

        for (int i=0; i<(int)receivedPackedMessages.size(); ++i){
            unpackMessages(receivedPackedMessages[i]);
            delete receivedPackedMessages[i];
        }
        receivedPackedMessages.clear();
    }

    combineMessages();
}

const char* Lumberjack::packMessages()
{
    //This function packs the messages into one long char array
    //  in the following format:
    //
    // <message count><largest message size><packed message size><packed message><packed message size>...
 
    if (m_messages.size() == 0) {
        return "0";
    }

    int totalSize = 1; // include size for null terminator

    //Calculate total size of char array after all messages are
    //  combined. Keep track largest message size so we don't have
    //  to realloc when unpacking
    std::vector<std::string> packedMessages;
    std::vector<std::string> sizeStrings;
    int currSize, largestSize = 0;
    int messageCount = (int)m_messages.size();
    for (int i=0;i<messageCount; ++i) {
        packedMessages.push_back(m_messages[i]->pack());
        currSize = packedMessages[i].size();
        sizeStrings.push_back(intToString(currSize));
        //           message size + size string size + memberDelimiter size
        totalSize += currSize + sizeStrings[i].size() + 1;
        if (largestSize < currSize) {
            largestSize = currSize;
        }
    }

    // Create and calculate size of message count and largest message size strings
    std::string largestSizeString = intToString(largestSize) + memberDelimiter;
    std::string messageCountString = intToString(messageCount) + memberDelimiter;
    totalSize += largestSizeString.size() + messageCountString.size();

    const char* packedMessagesString = new char[totalSize];
    char* packedMessagesIndex = (char*)packedMessagesString;

    // Copy largest size to start of packed message
    std::memcpy(packedMessagesIndex, messageCountString.c_str(), messageCountString.size());
    packedMessagesIndex += messageCountString.size();
    std::memcpy(packedMessagesIndex, largestSizeString.c_str(), largestSizeString.size());
    packedMessagesIndex += largestSizeString.size();

    for (int i=0;i<messageCount; ++i) {
        // Copy current message size
        std::memcpy(packedMessagesIndex, sizeStrings[i].c_str(), sizeStrings[i].size());
        packedMessagesIndex += sizeStrings[i].size();
        // Copy memberDelimiter
        std::memcpy(packedMessagesIndex, &memberDelimiter, sizeof(char));  //ToDo: better way to copy this im sure
        packedMessagesIndex += 1;
        // Copy packed message
        std::memcpy(packedMessagesIndex, packedMessages[i].c_str(), packedMessages[i].size());
        packedMessagesIndex += packedMessages[i].size();
    }
    packedMessagesIndex = '\0';
    return packedMessagesString;
}

void Lumberjack::unpackMessages(const char* packedMessages)
{
    int largestSize, messageCount;
    std::string tempString;
    int packedMessagesSize = std::strlen(packedMessages);
    int i = 0;
    // Get message count
    for(; i < packedMessagesSize; ++i) {
        if (packedMessages[i] == memberDelimiter) {
            messageCount = stringToInt(tempString);
            ++i;
            break;
        }
        tempString += packedMessages[i];
    }

    // Get largest message size
    tempString = "";
    for(; i < packedMessagesSize; ++i) {
        if (packedMessages[i] == memberDelimiter) {
            largestSize = stringToInt(tempString);
            ++i;
            break;
        }
        tempString += packedMessages[i];
    }

    // Grab each message    
    char* buffer = new char[largestSize+1];
    Message* message;
    int messageSize;
    for (int j = 0; j < messageCount; ++j) {
        tempString = "";
        for(; i < packedMessagesSize; ++i) {
            if (packedMessages[i] == memberDelimiter) {
                messageSize = stringToInt(tempString);
                ++i;
                break;
            }
            tempString += packedMessages[i];
        }
        memcpy(buffer, &packedMessages[i], messageSize*sizeof(char));
        buffer[messageSize] = '\0';
        message = new Message();
        message->unpack(std::string(buffer), m_ranksLimit);
        m_messages.push_back(message);
        i += messageSize;
    }
}

void Lumberjack::combineMessages()
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
