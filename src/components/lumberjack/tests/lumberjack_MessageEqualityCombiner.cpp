/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

#include "gtest/gtest.h"

#include "lumberjack/MessageEqualityCombiner.hpp"

TEST(lumberjack_MessageEqualityCombiner, case01)
{
    //Test most basic case: two equal messages
    std::string message = "I never wanted to do this job in the first place!";
    asctoolkit::lumberjack::MessageInfo m1;
    m1.message(message);
    m1.addRank(13, 5);
    m1.fileName("foo.cpp");
    m1.lineNumber(154);

    asctoolkit::lumberjack::MessageInfo m2;
    m2.message(message);
    m2.addRank(14, 5);
    m2.fileName("foo.cpp");
    m2.lineNumber(154);

    asctoolkit::lumberjack::MessageEqualityCombiner c;

    bool shouldMessageInfosBeCombined = c.shouldMessageInfosBeCombined(m1, m2);

    c.combine(m1, m2, 5);

    EXPECT_EQ(shouldMessageInfosBeCombined, true);
    EXPECT_EQ(m1.message().compare(message), 0);
    EXPECT_EQ(m1.rankCount(), 2);
    EXPECT_EQ(m1.ranks()[0], 13);
    EXPECT_EQ(m1.ranks()[1], 14);
}

TEST(lumberjack_MessageEqualityCombiner, case02)
{
    //Test most basic negative case: two not-equal messages
    std::string message1 = "I never wanted to do this job in the first place!";
    asctoolkit::lumberjack::MessageInfo m1;
    m1.message(message1);
    m1.addRank(13, 5);
    m1.fileName("foo.cpp");
    m1.lineNumber(154);

    std::string message2 = "This message is not equal to the first.";
    asctoolkit::lumberjack::MessageInfo m2;
    m2.message(message2);
    m2.addRank(14, 5);
    m2.fileName("foo.cpp");
    m2.lineNumber(154);

    asctoolkit::lumberjack::MessageEqualityCombiner c;

    bool shouldMessageInfosBeCombined = c.shouldMessageInfosBeCombined(m1, m2);

    EXPECT_EQ(shouldMessageInfosBeCombined, false);
    EXPECT_EQ(m1.message().compare(message1), 0);
    EXPECT_EQ(m1.rankCount(), 1);
    EXPECT_EQ(m1.ranks()[0], 13);

    EXPECT_EQ(m2.message().compare(message2), 0);
    EXPECT_EQ(m2.rankCount(), 1);
    EXPECT_EQ(m2.ranks()[0], 14);
}
