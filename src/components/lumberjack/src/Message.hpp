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
 * \file Message.hpp
 * \author Chris White (white238@llnl.gov)
 *
 * \brief This file contains the class definition of the Message.
 *******************************************************************************
 */

#ifndef MESSAGE_HPP
#define MESSAGE_HPP

#include <string>
#include <vector>

namespace asctoolkit {
namespace lumberjack {

/*!
 *******************************************************************************
 * \class Message
 *
 * \brief Holds all necessary information about messages and where they came from.
 *
 *  This class holds all the information about an individual message and where it
 *  came from, such as rank, file name, and line number.
 *
 * \see Combiner Logger
 *******************************************************************************
 */
class Message {
    public:
        // Constructors
          /*!
         *****************************************************************************
         * \brief Basic constructor where everything defaults to nothing.
         *****************************************************************************
         */
        Message()
        : m_text("")
        , m_ranksCount(0)
        , m_fileName("")
        , m_lineNumber(0) {}

        /*!
         *****************************************************************************
         * \brief Constructor where you can specify all values for a Message that originated
         * from a specific rank.
         *
         * \param [in] text Actual text of the Message.
         * \param [in] rank The rank where the Message originated.
         * \param [in] fileName The file name where the Message originated.
         * \param [in] lineNumber The line number where the Message originated.
         *****************************************************************************
         */
        Message(const std::string& text, int rank,
                    const std::string& fileName, int lineNumber)
        : m_text(text)
        , m_ranksCount(1)
        , m_fileName(fileName)
        , m_lineNumber(lineNumber)
        {
            m_ranks.push_back(rank);
        }

        /*!
         *****************************************************************************
         * \brief Constructor where you can specify all values for a Message that originated
         * from a multiple ranks.
         *
         * \param [in] text Actual text of the Message.
         * \param [in] ranks The rank where the Message originated.
         * \param [in] ranksCount Total amount of ranks where this Message has originated from.
         * \param [in] ranksLimit Limit on how many ranks are individually tracked per Message.
         * \param [in] fileName The file name where the Message originated.
         * \param [in] lineNumber The line number where the Message originated.
         *****************************************************************************
         */
        Message(const std::string& text, const std::vector<int>& ranks,
                    int ranksCount, int ranksLimit,
                    const std::string& fileName, int lineNumber)
        : m_text(text)
        , m_ranksCount(0)
        , m_fileName(fileName)
        , m_lineNumber(lineNumber)
        {
            addRanks(ranks, ranksCount, ranksLimit);
        }

        // Getters
        /*!
         *****************************************************************************
         * \brief Returns text of the Message.
         *****************************************************************************
         */
        std::string text() const;

        /*!
         *****************************************************************************
         * \brief Returns vector of the ranks where this Message originated.
         *****************************************************************************
         */
        std::vector<int> ranks() const;

        /*!
         *****************************************************************************
         * \brief Returns total count of ranks where this Message originated.
         *****************************************************************************
         */
        int ranksCount() const;

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
         * \brief Returns file name of where this Message originated.
         *****************************************************************************
         */
        std::string fileName() const;

        /*!
         *****************************************************************************
         * \brief Returns line number of where this Message originated.
         *****************************************************************************
         */
        int lineNumber() const;

        // Setters

        /*!
         *****************************************************************************
         * \brief Sets a new text for this Message.
         *
         * \param [in] newMessage The new text to be set for this Message.
         *****************************************************************************
         */
        void text(const std::string& newText);

        /*!
         *****************************************************************************
         * \brief Sets a new file name for this Message.
         *
         * \param [in] newFileName The new file name to be set for this Message.
         *****************************************************************************
         */
        void fileName(const std::string& newFileName);

        /*!
         *****************************************************************************
         * \brief Sets a new line number for this Message.
         *
         * \param [in] newLineNumber The delimiter used to separate the ranks in returned string.
         *****************************************************************************
         */
        void lineNumber(int newLineNumber);

        /*!
         *****************************************************************************
         * \brief Adds a rank to this Message.
         *
         * \param [in] newRank The new rank to be added.
         * \param [in] ranksLimit Limit on how many ranks are individually tracked per Message.
         *****************************************************************************
         */
        void addRank(int newRank, int ranksLimit);

        /*!
         *****************************************************************************
         * \brief Adds multiple ranks to this Message.  ranksCount is used to increment since
         *  duplicates are removed from Message::ranks.
         *
         * \param [in] newRanks The new ranks to be added.
         * \param [in] ranksCount Count to add to Message::ranksCount
         * \param [in] ranksLimit Limit on how many ranks are individually tracked per Message.
         *****************************************************************************
         */
        void addRanks(const std::vector<int>& newRanks, int ranksCount, int ranksLimit);

        // utilities

        /*!
         *****************************************************************************
         * \brief Returns a string of all information about this Message packed into a string.
         *
         * The Message is packed into a string utilizing the following format:
         *  <ranks delimited by ,>*<rank count>*<file name>*<line number>*<text>
         *
         *****************************************************************************
         */
        std::string pack();

        /*!
         *****************************************************************************
         * \brief Overrides the information in this Message with the given packed string.
         *
         * The Message is unpacked from a string utilizing the following format:
         *  <ranks delimited by ,>*<rank count>*<file name>*<line number>*<text>
         *
         * \param [in] packedMessage Packed Message containing the new information.
         * \param [in] ranksLimit The delimiter used to separate the ranks in returned string.
         *****************************************************************************
         */
        void unpack(const std::string& packedMessage, int ranksLimit);
    private:
        std::string m_text;
        std::vector<int> m_ranks;
        int m_ranksCount;
        std::string m_fileName;
        int m_lineNumber;
};

} // end namespace lumberjack
} // end namespace asctoolkit

#endif
