// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/slic/streams/LumberjackStream.hpp"

#include <vector>

#include "axom/core/Macros.hpp"

#include "axom/lumberjack/BinaryTreeCommunicator.hpp"
#include "axom/lumberjack/Lumberjack.hpp"

namespace axom
{
namespace slic
{
//------------------------------------------------------------------------------
LumberjackStream::LumberjackStream(std::ostream* stream,
                                   MPI_Comm comm,
                                   int ranksLimit)
  : m_isLJOwnedBySLIC(false)
  , m_stream(stream)
{
  this->initializeLumberjack(comm, ranksLimit);
}

//------------------------------------------------------------------------------
LumberjackStream::LumberjackStream(std::ostream* stream,
                                   MPI_Comm comm,
                                   int ranksLimit,
                                   const std::string& format)
  : m_isLJOwnedBySLIC(false)
  , m_stream(stream)
{
  this->initializeLumberjack(comm, ranksLimit);
  this->setFormatString(format);
}

//------------------------------------------------------------------------------
LumberjackStream::LumberjackStream(std::ostream* stream,
                                   axom::lumberjack::Lumberjack* lj)
  : m_lj(lj)
  , m_isLJOwnedBySLIC(false)
  , m_stream(stream)
{ }

//------------------------------------------------------------------------------
LumberjackStream::LumberjackStream(std::ostream* stream,
                                   axom::lumberjack::Lumberjack* lj,
                                   const std::string& format)
  : m_lj(lj)
  , m_isLJOwnedBySLIC(false)
  , m_stream(stream)
{
  this->setFormatString(format);
}

//------------------------------------------------------------------------------
LumberjackStream::~LumberjackStream()
{
  if(m_isLJOwnedBySLIC)
  {
    this->finalizeLumberjack();
  }
}

//------------------------------------------------------------------------------
void LumberjackStream::append(message::Level msgLevel,
                              const std::string& message,
                              const std::string& tagName,
                              const std::string& fileName,
                              int line,
                              bool AXOM_NOT_USED(filter_duplicates))
{
  if(m_lj == nullptr)
  {
    std::cerr
      << "ERROR: NULL Lumberjack instance in LumberjackStream::append!\n";
    return;
  }

  m_lj->queueMessage(message, fileName, line, msgLevel, tagName);
}

//------------------------------------------------------------------------------
void LumberjackStream::flush()
{
  if(m_lj == nullptr)
  {
    std::cerr
      << "ERROR: NULL Lumberjack instance in LumberjackStream::flush!\n";
    return;
  }

  m_lj->pushMessagesFully();
  this->write();
}

//------------------------------------------------------------------------------
void LumberjackStream::push()
{
  if(m_lj == nullptr)
  {
    std::cerr << "ERROR: NULL Lumberjack instance in LumberjackStream::push!\n";
    return;
  }

  m_lj->pushMessagesOnce();
}

//------------------------------------------------------------------------------
void LumberjackStream::write()
{
  if(m_lj == nullptr)
  {
    std::cerr
      << "ERROR: NULL Lumberjack instance in LumberjackStream::write!\n";
    return;
  }

  if(m_lj->isOutputNode())
  {
    std::vector<axom::lumberjack::Message*> messages = m_lj->getMessages();

    const int nmessages = static_cast<int>(messages.size());
    std::string rankString;
    for(int i = 0; i < nmessages; ++i)
    {
      rankString = std::to_string(messages[i]->count()) + ": " +
        messages[i]->stringOfRanks();
      (*m_stream) << this->getFormatedMessage(
        message::getLevelAsString(
          static_cast<message::Level>(messages[i]->level())),
        messages[i]->text(),
        messages[i]->tag(),
        rankString,
        messages[i]->fileName(),
        messages[i]->lineNumber());
    }

    m_lj->clearMessages();
  }
}

//------------------------------------------------------------------------------
void LumberjackStream::initializeLumberjack(MPI_Comm comm, int ranksLimit)
{
  m_ljComm = new axom::lumberjack::BinaryTreeCommunicator;
  m_ljComm->initialize(comm, ranksLimit);
  m_lj = new axom::lumberjack::Lumberjack;
  m_lj->initialize(m_ljComm, ranksLimit);
  m_isLJOwnedBySLIC = true;
}

//------------------------------------------------------------------------------
void LumberjackStream::finalizeLumberjack()
{
  m_lj->finalize();
  m_ljComm->finalize();
  delete m_lj;
  delete m_ljComm;
  m_isLJOwnedBySLIC = false;
}

} /* namespace slic */
} /* namespace axom */
