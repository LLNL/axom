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
 * \file Combiner.hpp
 * \author Chris White (white238@llnl.gov)
 *
 * \brief This file contains the abstract base class defining the interface of all Combiners.
*******************************************************************************
 */

#ifndef COMBINER_HPP
#define COMBINER_HPP

#include "lumberjack/MessageInfo.hpp"

namespace asctoolkit {
namespace lumberjack {

/*!
 *******************************************************************************
 * \class Combiner
 *
 * \brief Abstract base class defining the interface of all Combiner classes.
 *
 *  Concrete instances need to inherit from this class and implement these functions.
 *  You will need to add your Combiner using Logger::addCombiner
 *
 * \see MessageEqualityCombiner Logger
 *******************************************************************************
 */
class Combiner {
    public:
        /*!
         *****************************************************************************
         * \brief Returns the unique string identifier for this combiner. Used by Logger
         *  to differentiate between other combiners.
         *****************************************************************************
         */
        virtual const std::string id() = 0;

        /*!
         *****************************************************************************
         * \brief Function used by Logger to indicate whether two MessageInfo classes should be
         * combined.  They are not actually combined by this function.
         *
         * \param [in] leftMessageInfo The left MessageInfo to be compared.
         * \param [in] rightMessageInfo The right MessageInfo to be compared.
         *****************************************************************************
         */
        virtual bool shouldMessageInfosBeCombined(const MessageInfo& leftMessageInfo,
                                                  const MessageInfo& rightMessageInfo) = 0;

        /*!
         *****************************************************************************
         * \brief Combines the combinee into the combined MessageInfo.
         *
         * \param [in,out] combined the MessageInfo that will be modified.
         * \param [in] combinee the MessageInfo that is combined into the other.
         * \param [in] ranksLimit The limit on how many individual ranks are tracked in
         * the combined MessageInfo. MessageInfo::rankCount is always incremented.
         *****************************************************************************
         */
        virtual void combine(MessageInfo& combined,
                             const MessageInfo& combinee, const int ranksLimit) = 0;
};

} // end namespace lumberjack
} // end namespace asctoolkit

#endif
