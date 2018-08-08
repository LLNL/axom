/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*!
 ******************************************************************************
 *
 * \file BinaryTreeCommunicator.cpp
 *
 * \brief Implementation of the Lumberjack class.
 *
 ******************************************************************************
 */

#include "axom/lumberjack/Lumberjack.hpp"

#include "axom/core/utilities/StringUtilities.hpp"

#include <cstring>

namespace axom
{
namespace lumberjack
{

void Lumberjack::initialize(Communicator* communicator, int ranksLimit)
{
  m_communicator = communicator;
  m_ranksLimit = ranksLimit;
  m_combiners.push_back(new TextEqualityCombiner);
}

void Lumberjack::finalize()
{
  m_communicator = nullptr;
  clearCombiners();
  clearMessages();
}

void Lumberjack::addCombiner(Combiner* combiner)
{
  bool identifierFound = false;
  for (int i=0 ; i<(int)m_combiners.size() ; ++i)
  {
    if (m_combiners[i]->id() == combiner->id())
    {
      identifierFound = true;
      break;
    }
  }
  if (!identifierFound)
  {
    m_combiners.push_back(combiner);
  }
}

void Lumberjack::removeCombiner(const std::string& combinerIdentifier)
{
  int combinerToBeRemoved = -1;
  for (int i=0 ; i<(int)m_combiners.size() ; ++i)
  {
    if (m_combiners[i]->id() == combinerIdentifier)
    {
      delete m_combiners[i];
      combinerToBeRemoved = i;
      break;
    }
  }
  if (combinerToBeRemoved != -1)
  {
    m_combiners.erase(m_combiners.begin()+combinerToBeRemoved);
  }
}

void Lumberjack::clearCombiners()
{
  for (int i=0 ; i<(int)m_combiners.size() ; ++i)
  {
    delete m_combiners[i];
  }
  m_combiners.clear();
}

const std::vector<Message*>& Lumberjack::getMessages() const
{
  return m_messages;
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
  for (int i=0 ; i<(int)m_messages.size() ; ++i)
  {
    delete m_messages[i];
  }
  m_messages.clear();
}

void Lumberjack::queueMessage(const std::string& text)
{
  queueMessage(text, "", -1, 0, "");
}

void Lumberjack::queueMessage(const std::string& text,
                              const std::string& fileName, const int lineNumber,
                              int level, const std::string& tag)
{
  Message* mi = new Message(text, m_communicator->rank(),
                            fileName, lineNumber, level, tag);
  m_messages.push_back(mi);
}

void Lumberjack::pushMessagesOnce()
{
  const char* packedMessagesToBeSent = "";
  if (!m_communicator->isOutputNode())
  {
    combineMessages();
    packedMessagesToBeSent = packMessages();
    clearMessages();
  }
  std::vector<const char*> receivedPackedMessages;

  m_communicator->push(packedMessagesToBeSent, receivedPackedMessages);

  if (!m_communicator->isOutputNode())
  {
    delete [] packedMessagesToBeSent;
  }

  for (int i=0 ; i<(int)receivedPackedMessages.size() ; ++i)
  {
    unpackMessages(receivedPackedMessages[i]);
    delete [] receivedPackedMessages[i];
  }
  receivedPackedMessages.clear();

  combineMessages();
}

void Lumberjack::pushMessagesFully()
{
  const char* packedMessagesToBeSent = "";
  std::vector<const char*> receivedPackedMessages;
  int numPushesToFlush = m_communicator->numPushesToFlush();
  for (int i=0 ; i<numPushesToFlush ; ++i)
  {
    if (!m_communicator->isOutputNode())
    {
      combineMessages();
      packedMessagesToBeSent = packMessages();
      clearMessages();
    }

    m_communicator->push(packedMessagesToBeSent, receivedPackedMessages);

    if (!m_communicator->isOutputNode())
    {
      delete [] packedMessagesToBeSent;
    }

    for (int i=0 ; i<(int)receivedPackedMessages.size() ; ++i)
    {
      unpackMessages(receivedPackedMessages[i]);
      delete [] receivedPackedMessages[i];
    }
    receivedPackedMessages.clear();
  }

  combineMessages();
}

bool Lumberjack::isOutputNode()
{
  return m_communicator->isOutputNode();
}

const char* Lumberjack::packMessages()
{
  //This function packs the messages into one long char array
  //  in the following format:
  //
  // <message count><packed message size><packed message><packed message
  // size>...

  if (m_messages.size() == 0)
  {
    char* zero = new char[2];
    std::strcpy(zero, "0");
    return zero;
  }

  int totalSize = 1;   // include size for null terminator

  //Calculate total size of char array after all messages are
  //  combined. Keep track largest message size so we don't have
  //  to realloc when unpacking
  std::vector<std::string> packedMessages;
  std::vector<std::string> sizeStrings;
  int currSize, largestSize = 0;
  int messageCount = (int)m_messages.size();
  for (int i=0 ; i<messageCount ; ++i)
  {
    packedMessages.push_back(m_messages[i]->pack());
    currSize = packedMessages[i].size();
    sizeStrings.push_back(axom::utilities::string::intToString(currSize));
    //           message size + size string size + memberDelimiter size
    totalSize += currSize + sizeStrings[i].size() + 1;
    if (largestSize < currSize)
    {
      largestSize = currSize;
    }
  }

  // Create and calculate size of message count
  std::string messageCountString =
    axom::utilities::string::intToString(messageCount) + memberDelimiter;
  totalSize += messageCountString.size();

  const char* packedMessagesString = new char[totalSize];
  char* packedMessagesIndex = (char*)packedMessagesString;

  // Copy message count to start of packed message
  std::memcpy(packedMessagesIndex, messageCountString.c_str(),
              messageCountString.size());
  packedMessagesIndex += messageCountString.size();

  for (int i=0 ; i<messageCount ; ++i)
  {
    // Copy current message size
    std::memcpy(packedMessagesIndex, sizeStrings[i].c_str(),
                sizeStrings[i].size());
    packedMessagesIndex += sizeStrings[i].size();
    // Copy memberDelimiter
    std::memcpy(packedMessagesIndex, &memberDelimiter, sizeof(char));      //ToDo:
                                                                           // better
                                                                           // way
                                                                           // to
                                                                           // copy
                                                                           // this
                                                                           // im
                                                                           // sure
    packedMessagesIndex += 1;
    // Copy packed message
    std::memcpy(packedMessagesIndex, packedMessages[i].c_str(),
                packedMessages[i].size());
    packedMessagesIndex += packedMessages[i].size();
  }

  packedMessagesIndex[0] = '\0';
  return packedMessagesString;
}

void Lumberjack::unpackMessages(const char* packedMessages)
{
  std::string packedMessagesString = std::string(packedMessages);
  std::size_t start, end;
  std::string tempSubString = "";

  // Get message count
  end = packedMessagesString.find(memberDelimiter);
  tempSubString = packedMessagesString.substr(0, end);
  int messageCount = axom::utilities::string::stringToInt(tempSubString);
  start = end + 1;

  // Grab each message
  Message* message;
  int messageSize;
  for (int j = 0 ; j < messageCount ; ++j)
  {
    //Get current message size
    end = packedMessagesString.find(memberDelimiter, start);
    tempSubString = packedMessagesString.substr(start, end-start);
    messageSize = axom::utilities::string::stringToInt(tempSubString);
    start = end + 1;

    //Create current message and save
    message = new Message();
    tempSubString = packedMessagesString.substr(start, messageSize);
    message->unpack(tempSubString, m_ranksLimit);
    m_messages.push_back(message);
    start += messageSize;
  }
}

void Lumberjack::combineMessages()
{
  int messagesSize = (int)m_messages.size();
  if (messagesSize < 2)
  {
    return;
  }

  std::vector<Message*> finalMessages;
  std::vector<int> indexesToBeDeleted;
  int combinersSize = (int)m_combiners.size();
  bool combinedMessage = false;
  finalMessages.push_back(m_messages[0]);
  for (int allIndex=1 ; allIndex<messagesSize ; ++allIndex)
  {
    combinedMessage = false;
    for (int finalIndex=0 ; finalIndex<(int)finalMessages.size() ; ++finalIndex)
    {
      for (int combinerIndex=0 ; combinerIndex<combinersSize ; ++combinerIndex)
      {
        if (m_combiners[combinerIndex]->
            shouldMessagesBeCombined(*finalMessages[finalIndex],
                                     *m_messages[allIndex]))
        {
          m_combiners[combinerIndex]->combine(*finalMessages[finalIndex],
                                              *m_messages[allIndex],
                                              m_ranksLimit);
          indexesToBeDeleted.push_back(allIndex);
          combinedMessage = true;
          break;
        }
      }
      if (combinedMessage)
      {
        break;
      }
    }
    if (!combinedMessage)
    {
      finalMessages.push_back(m_messages[allIndex]);
    }
  }

  for (int i=0 ; i<(int)indexesToBeDeleted.size() ; ++i)
  {
    delete m_messages[indexesToBeDeleted[i]];
  }
  m_messages.swap(finalMessages);
}
} // end namespace lumberjack
} // end namespace axom
