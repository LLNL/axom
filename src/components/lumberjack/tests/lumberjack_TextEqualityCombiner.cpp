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

#include "lumberjack/TextEqualityCombiner.hpp"

TEST(lumberjack_TextEqualityCombiner, case01)
{
    //Test most basic case: two equal texts
    std::string text = "I never wanted to do this job in the first place!";
    asctoolkit::lumberjack::Message m1;
    m1.text(text);
    m1.addRank(13, 5);
    m1.fileName("foo.cpp");
    m1.lineNumber(154);

    asctoolkit::lumberjack::Message m2;
    m2.text(text);
    m2.addRank(14, 5);
    m2.fileName("foo.cpp");
    m2.lineNumber(154);

    asctoolkit::lumberjack::TextEqualityCombiner c;

    bool shouldMessagesBeCombined = c.shouldMessagesBeCombined(m1, m2);

    c.combine(m1, m2, 5);

    EXPECT_EQ(shouldMessagesBeCombined, true);
    EXPECT_EQ(m1.text().compare(text), 0);
    EXPECT_EQ(m1.ranksCount(), 2);
    EXPECT_EQ(m1.ranks()[0], 13);
    EXPECT_EQ(m1.ranks()[1], 14);
}

TEST(lumberjack_TextEqualityCombiner, case02)
{
    //Test most basic negative case: two not-equal texts
    std::string text1 = "I never wanted to do this job in the first place!";
    asctoolkit::lumberjack::Message m1;
    m1.text(text1);
    m1.addRank(13, 5);
    m1.fileName("foo.cpp");
    m1.lineNumber(154);

    std::string text2 = "This text is not equal to the first.";
    asctoolkit::lumberjack::Message m2;
    m2.text(text2);
    m2.addRank(14, 5);
    m2.fileName("foo.cpp");
    m2.lineNumber(154);

    asctoolkit::lumberjack::TextEqualityCombiner c;

    bool shouldMessagesBeCombined = c.shouldMessagesBeCombined(m1, m2);

    EXPECT_EQ(shouldMessagesBeCombined, false);
    EXPECT_EQ(m1.text().compare(text1), 0);
    EXPECT_EQ(m1.ranksCount(), 1);
    EXPECT_EQ(m1.ranks()[0], 13);

    EXPECT_EQ(m2.text().compare(text2), 0);
    EXPECT_EQ(m2.ranksCount(), 1);
    EXPECT_EQ(m2.ranks()[0], 14);
}
