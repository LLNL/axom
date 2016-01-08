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
 * \file Message.cpp
 * \author Chris White (white238@llnl.gov)
 *
 * \brief This file contains the class implementation of the Message.
 *******************************************************************************
 */

#include "lumberjack/Message.hpp"

#include "lumberjack/Utility.hpp"

#include <algorithm>
#include <iostream>

namespace asctoolkit {
namespace lumberjack {

std::string Message::text() const
{
    return m_text;
}

void Message::text(const std::string& newText)
{
    m_text = newText;
}

std::vector<int> Message::ranks() const
{
    return m_ranks;
}

int Message::ranksCount() const
{
    return m_ranksCount;
}

std::string Message::fileName() const
{
    return m_fileName;
}

std::string Message::stringOfRanks(std::string delimiter) const
{
    std::string returnString = "";
    int ranksSize = m_ranks.size();
    for(int i=0; i<ranksSize;++i){
        returnString += intToString(m_ranks[i]);
        if (i < (ranksSize-1)) {
            returnString += delimiter;
        }
    }
    return returnString;
}

void Message::fileName(const std::string& newFileName)
{
    m_fileName = newFileName;
}

int Message::lineNumber() const
{
    return m_lineNumber;
}

void Message::lineNumber(int newLineNumber)
{
    m_lineNumber = newLineNumber;
}

void Message::addRank(int newRank, int ranksLimit)
{
    // If ranksLimit has already been reached don't add newRank to m_ranks
    if (m_ranks.size() < (std::vector<int>::size_type)ranksLimit){
        // If newRank is already in m_ranks then don't add it
        std::vector<int>::iterator iter = std::find(m_ranks.begin(), m_ranks.end(), newRank);
        if ((m_ranks.size() == 0) || (iter == m_ranks.end())){
            m_ranks.push_back(newRank);
        }
    }
    // Always increment rank count
    m_ranksCount++;
}

void Message::addRanks(const std::vector<int>& newRanks, int ranksCount, int ranksLimit)
{
    int newRanksSize = newRanks.size();
    for(int i=0; i<newRanksSize; ++i){
        // If ranksLimit has already been reached don't add newRank to m_ranks
        if (m_ranks.size() >= (std::vector<int>::size_type)ranksLimit){
            break;
        }
        // If newRank is already in m_ranks then don't add it
        std::vector<int>::iterator iter = std::find(m_ranks.begin(), m_ranks.end(), newRanks[i]);
        if ((m_ranks.size() == 0) || (iter == m_ranks.end())){
            m_ranks.push_back(newRanks[i]);
        }
    }
    // Always increment ranks count
    m_ranksCount += ranksCount;
}

std::string Message::pack()
{
    std::string packedMessage;
    int ranksSize = (int)m_ranks.size();
    for (int i=0; i<ranksSize; ++i){
        packedMessage += intToString(m_ranks[i]);
        if (i < (ranksSize-1)) {
            packedMessage += rankDelimiter;
        }
    }
    packedMessage += memberDelimiter + intToString(m_ranksCount);
    packedMessage += memberDelimiter + m_fileName + memberDelimiter;

    if (m_lineNumber > 0){
        packedMessage += intToString(m_lineNumber);
    }
    packedMessage += memberDelimiter + m_text;

    return packedMessage;
}

void Message::unpack(const std::string& packedMessage, int ranksLimit)
{
    int messageLength = (int)packedMessage.length();
    std::string currString;
    int i = 0;

    // Grab ranks
    m_ranks.clear();
    int currRank;
    for (i=0; i<messageLength; ++i) {
        if ((packedMessage[i] == memberDelimiter) ||
            (packedMessage[i] == rankDelimiter)) {
            currRank = stringToInt(currString);
            currString = "";
            addRank(currRank, ranksLimit);
            if (packedMessage[i] == memberDelimiter) {
                ++i;
                break;
            }
            else {
                continue;
            }
        }
        currString += packedMessage[i];
    }
    if (i >= messageLength) {
        //ToDo: figure out a better error handling method
        std::cout << "Error: Lumberjack recieved a truncated message that ended in the rank section." << std::endl;
        std::cout << packedMessage << std::endl;
    }

    //Grab rank count since it can differ from list that is sent
    currString = "";
    for (; i<messageLength; ++i) {
        if (packedMessage[i] == memberDelimiter) {
            m_ranksCount = stringToInt(currString);
            ++i;
            break;
        }
        currString += packedMessage[i];
    }
    if (i >= messageLength) {
        //ToDo: figure out a better error handling method
        std::cout << "Error: Lumberjack recieved a truncated message that ended in the rank count section." << std::endl;
        std::cout << packedMessage << std::endl;
    }

    //Grab file name
    currString = "";
    for (; i<messageLength; ++i) {
        if (packedMessage[i] == memberDelimiter) {
            m_fileName = currString;
            ++i;
            break;
        }
        currString += packedMessage[i];
    }
    if (i >= messageLength) {
        //ToDo: figure out a better error handling method
        std::cout << "Error: Lumberjack recieved a truncated message that ended in the file name section." << std::endl;
        std::cout << packedMessage << std::endl;
    }

    //Grab line number
    currString = "";
    for (; i<messageLength; ++i) {
        if (packedMessage[i] == memberDelimiter) {
            m_lineNumber = stringToInt(currString);
            ++i;
            break;
        }
        currString += packedMessage[i];
    }
    if (i >= messageLength) {
        //ToDo: figure out a better error handling method
        std::cout << "Error: Lumberjack recieved a truncated message that ended in the line number section." << std::endl;
        std::cout << packedMessage << std::endl;
    }

    //Grab message
    currString = "";
    for (; i<messageLength; ++i) {
        currString += packedMessage[i];
    }
    m_text = currString;
}

} // end namespace lumberjack
} // end namespace asctoolkit
