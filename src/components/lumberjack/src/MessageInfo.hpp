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
 * \file MessageInfo.hpp
 * \author Chris White (white238@llnl.gov)
 *
 * \brief This file contains the class definition of the MessageInfo.
 *******************************************************************************
 */

#ifndef MESSAGEINFO_HPP
#define MESSAGEINFO_HPP

#include <string>
#include <vector>

namespace asctoolkit {
namespace lumberjack {

/*!
 *******************************************************************************
 * \class MessageInfo
 *
 * \brief Holds all necessary information about messages and where they came from.
 *
 *  This class holds all the information about an individual message and where it
 *  came from, such as rank, file name, and line number.
 *
 * \see Combiner Logger
 *******************************************************************************
 */
class MessageInfo {
    public:
        // Constructors
          /*!
         *****************************************************************************
         * \brief Basic constructor where everything defaults to nothing.
         *****************************************************************************
         */
        MessageInfo()
        : m_message("")
        , m_rankCount(0)
        , m_fileName("")
        , m_lineNumber(0) {}

        /*!
         *****************************************************************************
         * \brief Constructor where you can specify all values for a message that originated
         * from a specific rank.
         *
         * \param [in] message Actual text of the message.
         * \param [in] rank The rank where the message originated.
         * \param [in] fileName The file name where the message originated.
         * \param [in] lineNumber The line number where the message originated.
         *****************************************************************************
         */
        MessageInfo(const std::string& message, int rank,
                    const std::string& fileName, int lineNumber)
        : m_message(message)
        , m_rankCount(1)
        , m_fileName(fileName)
        , m_lineNumber(lineNumber)
        {
            m_ranks.push_back(rank);
        }

        /*!
         *****************************************************************************
         * \brief Constructor where you can specify all values for a message that originated
         * from a multiple ranks.
         *
         * \param [in] message Actual text of the message.
         * \param [in] ranks The rank where the message originated.
         * \param [in] rankCount Total amount of ranks where this message has originated from.
         * \param [in] ranksLimit Limit on how many ranks are individually tracked per MessageInfo.
         * \param [in] fileName The file name where the message originated.
         * \param [in] lineNumber The line number where the message originated.
         *****************************************************************************
         */
        MessageInfo(const std::string& message, const std::vector<int>& ranks,
                    int rankCount, int ranksLimit,
                    const std::string& fileName, int lineNumber)
        : m_message(message)
        , m_fileName(fileName)
        , m_lineNumber(lineNumber)
        {
            addRanks(ranks, ranksLimit);
            m_rankCount = rankCount;
        }

        // Getters
        /*!
         *****************************************************************************
         * \brief Returns text of the message.
         *****************************************************************************
         */
        std::string message() const;

        /*!
         *****************************************************************************
         * \brief Returns vector of the ranks where this message originated.
         *****************************************************************************
         */
        std::vector<int> ranks() const;

        /*!
         *****************************************************************************
         * \brief Returns total rank count of where this message originated.
         *****************************************************************************
         */
        int rankCount() const;

        /*!
         *****************************************************************************
         * \brief Returns a string of ranks delimited by ',' unless otherwise specified.
         *
         * \param [in] delimiter The delimiter used to separate the ranks in returned string.
         *****************************************************************************
         */
        std::string stringOfRanks(std::string delimiter=",") const;

        /*!
         *****************************************************************************
         * \brief Returns file name of where this message originated.
         *****************************************************************************
         */
        std::string fileName() const;

        /*!
         *****************************************************************************
         * \brief Returns line number of where this message originated.
         *****************************************************************************
         */
        int lineNumber() const;

        // Setters

        /*!
         *****************************************************************************
         * \brief Sets a new message for this MessageInfo.
         *
         * \param [in] newMessage The new message to be set for this message.
         *****************************************************************************
         */
        void message(const std::string& newMessage);

        /*!
         *****************************************************************************
         * \brief Sets a new file name for this MessageInfo.
         *
         * \param [in] newFileName The new file name to be set for this message.
         *****************************************************************************
         */
        void fileName(const std::string& newFileName);

        /*!
         *****************************************************************************
         * \brief Sets a new line number for this MessageInfo.
         *
         * \param [in] newLineNumber The delimiter used to separate the ranks in returned string.
         *****************************************************************************
         */
        void lineNumber(int newLineNumber);

        /*!
         *****************************************************************************
         * \brief Adds a rank to this MessageInfo.
         *
         * \param [in] newRank The new rank to be added.
         * \param [in] ranksLimit Limit on how many ranks are individually tracked per MessageInfo.
         *****************************************************************************
         */
        void addRank(int newRank, int ranksLimit);

        /*!
         *****************************************************************************
         * \brief Adds multiple ranks to this MessageInfo.
         *
         * \param [in] newRanks The new ranks to be added.
         * \param [in] ranksLimit Limit on how many ranks are individually tracked per MessageInfo.
         *****************************************************************************
         */
        void addRanks(const std::vector<int>& newRanks, int ranksLimit);

        // utilities

        /*!
         *****************************************************************************
         * \brief Returns a string of all information about this message packed into a string.
         *
         * The MessageInfo is packed into a string utilizing the following format:
         *  <ranks delimited by ,>*<rank count>*<file name>*<line number>*<message>
         *
         *****************************************************************************
         */
        std::string pack();

        /*!
         *****************************************************************************
         * \brief Overrides the information in this MessageInfo with the given packed string.
         *
         * The MessageInfo is unpacked from a string utilizing the following format:
         *  <ranks delimited by ,>*<rank count>*<file name>*<line number>*<message>
         *
         * \param [in] packedMessage Packed message containing the new information.
         * \param [in] ranksLimit The delimiter used to separate the ranks in returned string.
         *****************************************************************************
         */
        void unpack(const std::string& packedMessage, int ranksLimit);
    private:
        std::string m_message;
        std::vector<int> m_ranks;
        int m_rankCount;
        std::string m_fileName;
        int m_lineNumber;
};

} // end namespace lumberjack
} // end namespace asctoolkit

#endif
