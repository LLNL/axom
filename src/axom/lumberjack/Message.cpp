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
 * \file Message.cpp
 *
 * \brief Implementation of the Message class.
 *
 ******************************************************************************
 */

#include "axom/lumberjack/Message.hpp"

#include "axom/core/utilities/StringUtilities.hpp"

#include <algorithm>
#include <iostream>

namespace axom
{
namespace lumberjack
{

//Getters

std::string Message::text() const
{
  return m_text;
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

int Message::lineNumber() const
{
  return m_lineNumber;
}

int Message::level() const
{
  return m_level;
}

std::string Message::tag() const
{
  return m_tag;
}

std::string Message::stringOfRanks(std::string delimiter) const
{
  std::string returnString = "";
  int ranksSize = m_ranks.size();
  for(int i=0 ; i<ranksSize ; ++i)
  {
    returnString += axom::utilities::string::intToString(m_ranks[i]);
    if (i < (ranksSize-1))
    {
      returnString += delimiter;
    }
  }
  return returnString;
}

//Setters

void Message::text(const std::string& newText)
{
  m_text = newText;
}

void Message::fileName(const std::string& newFileName)
{
  m_fileName = newFileName;
}

void Message::lineNumber(int newLineNumber)
{
  m_lineNumber = newLineNumber;
}

void Message::level(int newLevel)
{
  m_level = newLevel;
}

void Message::tag(const std::string& newTag)
{
  m_tag = newTag;
}

void Message::addRank(int newRank, int ranksLimit)
{
  // If ranksLimit has already been reached don't add newRank to m_ranks
  if (m_ranks.size() < (std::vector<int>::size_type)ranksLimit)
  {
    // If newRank is already in m_ranks then don't add it
    std::vector<int>::iterator iter = std::find(m_ranks.begin(),
                                                m_ranks.end(), newRank);
    if ((m_ranks.size() == 0) || (iter == m_ranks.end()))
    {
      m_ranks.push_back(newRank);
    }
  }
  // Always increment rank count
  m_ranksCount++;
}

void Message::addRanks(const std::vector<int>& newRanks, int ranksCount,
                       int ranksLimit)
{
  int newRanksSize = newRanks.size();
  for(int i=0 ; i<newRanksSize ; ++i)
  {
    // If ranksLimit has already been reached don't add newRank to m_ranks
    if (m_ranks.size() >= (std::vector<int>::size_type)ranksLimit)
    {
      break;
    }
    // If newRank is already in m_ranks then don't add it
    std::vector<int>::iterator iter = std::find(m_ranks.begin(),
                                                m_ranks.end(), newRanks[i]);
    if ((m_ranks.size() == 0) || (iter == m_ranks.end()))
    {
      m_ranks.push_back(newRanks[i]);
    }
  }
  // Always increment ranks count
  m_ranksCount += ranksCount;
}

// Utilities

std::string Message::pack()
{
  std::string packedMessage;

  int ranksSize = (int)m_ranks.size();
  for (int i=0 ; i<ranksSize ; ++i)
  {
    packedMessage += axom::utilities::string::intToString(m_ranks[i]);
    if (i < (ranksSize-1))
    {
      packedMessage += rankDelimiter;
    }
  }
  packedMessage += memberDelimiter;

  packedMessage += axom::utilities::string::intToString(m_ranksCount) +
                   memberDelimiter;

  packedMessage += m_fileName + memberDelimiter;

  if (m_lineNumber > 0)
  {
    packedMessage += axom::utilities::string::intToString(m_lineNumber);
  }
  packedMessage += memberDelimiter;

  packedMessage += axom::utilities::string::intToString(m_level) +
                   memberDelimiter;

  packedMessage += m_tag + memberDelimiter;

  packedMessage += m_text;

  return packedMessage;
}

void Message::unpack(const std::string& packedMessage, int ranksLimit)
{
  std::size_t start, end;
  std::string section;

  // Grab ranks
  end = packedMessage.find(memberDelimiter);
  if (end == std::string::npos)
  {
    std::cerr << "Error: Lumberjack recieved a truncated message "
              << "that ended in the ranks section."
              << std::endl;
    std::cerr << packedMessage << std::endl;
  }
  unpackRanks(packedMessage.substr(0, end), ranksLimit);
  start = end + 1;

  //Grab rank count since it can differ from list that is sent
  end = packedMessage.find(memberDelimiter, start);
  if (end == std::string::npos)
  {
    std::cerr << "Error: Lumberjack recieved a truncated message "
              << "that ended in the rank count section."
              << std::endl;
    std::cerr << packedMessage << std::endl;
  }
  m_ranksCount =
    axom::utilities::string::stringToInt(packedMessage.substr(start,
                                                              end-start));
  start = end + 1;

  //Grab file name
  end = packedMessage.find(memberDelimiter, start);
  if (end == std::string::npos)
  {
    std::cerr << "Error: Lumberjack recieved a truncated message "
              << "that ended in the file name section."
              << std::endl;
    std::cerr << packedMessage << std::endl;
  }
  m_fileName = packedMessage.substr(start, end-start);
  start = end + 1;

  //Grab line number
  end = packedMessage.find(memberDelimiter, start);
  if (end == std::string::npos)
  {
    std::cerr << "Error: Lumberjack recieved a truncated message "
              << "that ended in the line number section."
              << std::endl;
    std::cerr << packedMessage << std::endl;
  }
  if (end-start > 0)
  {
    m_lineNumber =
      axom::utilities::string::stringToInt(packedMessage.substr(start,
                                                                end-start));
  }
  start = end + 1;

  //Grab level
  end = packedMessage.find(memberDelimiter, start);
  if (end == std::string::npos)
  {
    std::cerr << "Error: Lumberjack recieved a truncated message "
              << "that ended in the level section."
              << std::endl;
    std::cerr << packedMessage << std::endl;
  }
  m_level =
    axom::utilities::string::stringToInt(packedMessage.substr(start,
                                                              end-start));
  start = end + 1;

  //Grab tag
  end = packedMessage.find(memberDelimiter, start);
  if (end == std::string::npos)
  {
    std::cerr << "Error: Lumberjack recieved a truncated message "
              << "that ended in the tag section."
              << std::endl;
    std::cerr << packedMessage << std::endl;
  }
  m_tag = packedMessage.substr(start, end-start);
  start = end + 1;

  //Grab message
  m_text = packedMessage.substr(end+1);
}

void Message::unpackRanks(const std::string& ranksString, int ranksLimit)
{
  m_ranks.clear();
  if (ranksString.empty())
  {
    std::cerr << "Error: Lumberjack recieved an empty rank section."
              << std::endl;
    return;
  }

  std::size_t start, end = ranksString.find(rankDelimiter);
  start = 0;
  while (true)
  {
    if (end == std::string::npos)
    {
      addRank(axom::utilities::string::stringToInt(ranksString.substr(
                                                     start)), ranksLimit);
      break;
    }
    else
    {
      addRank(axom::utilities::string::stringToInt(ranksString.substr(start,
                                                                      end-start)),
              ranksLimit);
    }
    start = end + 1;
    end = ranksString.find(rankDelimiter, start);
  }
}

} // end namespace lumberjack
} // end namespace axom
