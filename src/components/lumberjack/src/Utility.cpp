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
 * \file Utility.cpp
 * \author Chris White (white238@llnl.gov)
 *
 * \brief This file contains the implementation of utility functions.
 *******************************************************************************
 */

#include "lumberjack/Utility.hpp"

#ifndef ENABLE_CXX11
#include <sstream>
#endif

namespace asctoolkit {
namespace lumberjack {

void combineMessages(std::vector<Message*>& messages, std::vector<Combiner*>& combiners, int ranksLimit)
{
    int messagesSize = (int)messages.size();
    if (messagesSize < 2){
        return;
    }

    std::vector<Message*> finalMessages;
    std::vector<int> indexesToBeDeleted;
    int combinersSize = (int)combiners.size();
    bool combinedMessage = false;
    finalMessages.push_back(messages[0]);
    for (int allIndex=1; allIndex<messagesSize; ++allIndex){
        combinedMessage = false;
        for (int finalIndex=0; finalIndex<(int)finalMessages.size(); ++finalIndex){
            for (int combinerIndex=0; combinerIndex<combinersSize; ++combinerIndex){
                if (combiners[combinerIndex]->shouldMessagesBeCombined(*finalMessages[finalIndex],
                                                                       *messages[allIndex])){
                    combiners[combinerIndex]->combine(*finalMessages[finalIndex],
                                                      *messages[allIndex], ranksLimit);
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
            finalMessages.push_back(messages[allIndex]);
        }
    }

    for (int i=0; i<(int)indexesToBeDeleted.size(); ++i){
        delete messages[indexesToBeDeleted[i]];
    }
    messages.swap(finalMessages);
}

std::string intToString(int intValue)
{
    std::string stringValue = "";
#ifdef ENABLE_CXX11
    stringValue += std::to_string(intValue);
#else
    std::ostringstream ss;
    ss << intValue;
    stringValue += ss.str();
#endif
    return stringValue;
}

int stringToInt(const std::string& stringValue)
{
    int intValue = 0;
#ifdef ENABLE_CXX11
    intValue = stoi(stringValue);
#else
    std::istringstream(stringValue) >> intValue;
#endif
    return intValue;
}

} // end namespace lumberjack
} // end namespace asctoolkit
