// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/lumberjack/LineFileTagCombiner.hpp"

TEST(lumberjack_LineFileTagCombiner, case01)
{
  //Test positive case: two equal line numbers, filenames, and tags
  std::string text1 =
    "This message does not matter because we do not filter by text";
  axom::lumberjack::Message m1;
  m1.text(text1);
  m1.addRank(13, 5);
  m1.fileName("foo.cpp");
  m1.lineNumber(154);
  m1.tag("myTag");

  std::string text2 =
    "This message ALSO does not matter because we do not filter by text";
  axom::lumberjack::Message m2;
  m2.text(text2);
  m2.addRank(14, 5);
  m2.fileName("foo.cpp");
  m2.lineNumber(154);
  m2.tag("myTag");

  axom::lumberjack::LineFileTagCombiner c;

  bool shouldMessagesBeCombined = c.shouldMessagesBeCombined(m1, m2);

  if(shouldMessagesBeCombined == true)
  {
    c.combine(m1, m2, 5);
  }

  EXPECT_EQ(shouldMessagesBeCombined, true);

  EXPECT_EQ(m1.lineNumber(), 154);
  EXPECT_EQ(m1.fileName().compare("foo.cpp"), 0);
  EXPECT_EQ(m1.text().compare(text1), 0);
  EXPECT_EQ(m1.count(), 2);
  EXPECT_EQ(m1.ranks().size(), 2);
  EXPECT_EQ(m1.ranks()[0], 13);
  EXPECT_EQ(m1.ranks()[1], 14);
  EXPECT_EQ(m1.tag().compare("myTag"), 0);

  EXPECT_EQ(m2.lineNumber(), 154);
  EXPECT_EQ(m2.fileName().compare("foo.cpp"), 0);
  EXPECT_EQ(m2.text().compare(text2), 0);
  EXPECT_EQ(m2.count(), 1);
  EXPECT_EQ(m2.ranks().size(), 1);
  EXPECT_EQ(m2.ranks()[0], 14);
  EXPECT_EQ(m2.tag().compare("myTag"), 0);
}

TEST(lumberjack_LineFileTagCombiner, case02)
{
  //Test negative case: two not equal line numbers, equal filenames, and equal tags
  std::string text1 =
    "This message does not matter because we do not filter by text";
  axom::lumberjack::Message m1;
  m1.text(text1);
  m1.addRank(13, 5);
  m1.fileName("foo.cpp");
  m1.lineNumber(12000);
  m1.tag("myTag");

  std::string text2 =
    "This message ALSO does not matter because we do not filter by text";
  axom::lumberjack::Message m2;
  m2.text(text2);
  m2.addRank(14, 5);
  m2.fileName("foo.cpp");
  m2.lineNumber(154);
  m2.tag("myTag");

  axom::lumberjack::LineFileTagCombiner c;

  bool shouldMessagesBeCombined = c.shouldMessagesBeCombined(m1, m2);

  if(shouldMessagesBeCombined == true)
  {
    c.combine(m1, m2, 5);
  }

  EXPECT_EQ(shouldMessagesBeCombined, false);

  EXPECT_EQ(m1.lineNumber(), 12000);
  EXPECT_EQ(m1.fileName().compare("foo.cpp"), 0);
  EXPECT_EQ(m1.text().compare(text1), 0);
  EXPECT_EQ(m1.count(), 1);
  EXPECT_EQ(m1.ranks().size(), 1);
  EXPECT_EQ(m1.ranks()[0], 13);
  EXPECT_EQ(m1.tag().compare("myTag"), 0);

  EXPECT_EQ(m2.lineNumber(), 154);
  EXPECT_EQ(m2.fileName().compare("foo.cpp"), 0);
  EXPECT_EQ(m2.text().compare(text2), 0);
  EXPECT_EQ(m2.count(), 1);
  EXPECT_EQ(m2.ranks().size(), 1);
  EXPECT_EQ(m2.ranks()[0], 14);
  EXPECT_EQ(m2.tag().compare("myTag"), 0);
}

TEST(lumberjack_LineFileTagCombiner, case03)
{
  //Test negative case: two equal line numbers, not equal filenames, and equal tags
  std::string text1 =
    "This message does not matter because we do not filter by text";
  axom::lumberjack::Message m1;
  m1.text(text1);
  m1.addRank(13, 5);
  m1.fileName("foo.cpp");
  m1.lineNumber(154);
  m1.tag("myTag");

  std::string text2 =
    "This message ALSO does not matter because we do not filter by text";
  axom::lumberjack::Message m2;
  m2.text(text2);
  m2.addRank(14, 5);
  m2.fileName("bar.cpp");
  m2.lineNumber(154);
  m2.tag("myTag");

  axom::lumberjack::LineFileTagCombiner c;

  bool shouldMessagesBeCombined = c.shouldMessagesBeCombined(m1, m2);

  if(shouldMessagesBeCombined == true)
  {
    c.combine(m1, m2, 5);
  }

  EXPECT_EQ(shouldMessagesBeCombined, false);

  EXPECT_EQ(m1.lineNumber(), 154);
  EXPECT_EQ(m1.fileName().compare("foo.cpp"), 0);
  EXPECT_EQ(m1.text().compare(text1), 0);
  EXPECT_EQ(m1.count(), 1);
  EXPECT_EQ(m1.ranks().size(), 1);
  EXPECT_EQ(m1.ranks()[0], 13);
  EXPECT_EQ(m1.tag().compare("myTag"), 0);

  EXPECT_EQ(m2.lineNumber(), 154);
  EXPECT_EQ(m2.fileName().compare("bar.cpp"), 0);
  EXPECT_EQ(m2.text().compare(text2), 0);
  EXPECT_EQ(m2.count(), 1);
  EXPECT_EQ(m2.ranks().size(), 1);
  EXPECT_EQ(m2.ranks()[0], 14);
  EXPECT_EQ(m2.tag().compare("myTag"), 0);
}

TEST(lumberjack_LineFileTagCombiner, case04)
{
  //Test negative case: two equal line numbers, equal filenames, and not equal tags
  std::string text1 =
    "This message does not matter because we do not filter by text";
  axom::lumberjack::Message m1;
  m1.text(text1);
  m1.addRank(13, 5);
  m1.fileName("foo.cpp");
  m1.lineNumber(154);
  m1.tag("myTag");

  std::string text2 =
    "This message ALSO does not matter because we do not filter by text";
  axom::lumberjack::Message m2;
  m2.text(text2);
  m2.addRank(14, 5);
  m2.fileName("foo.cpp");
  m2.lineNumber(154);
  m2.tag("myOtherTag");

  axom::lumberjack::LineFileTagCombiner c;

  bool shouldMessagesBeCombined = c.shouldMessagesBeCombined(m1, m2);

  if(shouldMessagesBeCombined == true)
  {
    c.combine(m1, m2, 5);
  }

  EXPECT_EQ(shouldMessagesBeCombined, false);

  EXPECT_EQ(m1.lineNumber(), 154);
  EXPECT_EQ(m1.fileName().compare("foo.cpp"), 0);
  EXPECT_EQ(m1.text().compare(text1), 0);
  EXPECT_EQ(m1.count(), 1);
  EXPECT_EQ(m1.ranks().size(), 1);
  EXPECT_EQ(m1.ranks()[0], 13);
  EXPECT_EQ(m1.tag().compare("myTag"), 0);

  EXPECT_EQ(m2.lineNumber(), 154);
  EXPECT_EQ(m2.fileName().compare("foo.cpp"), 0);
  EXPECT_EQ(m2.text().compare(text2), 0);
  EXPECT_EQ(m2.count(), 1);
  EXPECT_EQ(m2.ranks().size(), 1);
  EXPECT_EQ(m2.ranks()[0], 14);
  EXPECT_EQ(m2.tag().compare("myOtherTag"), 0);
}

TEST(lumberjack_LineFileTagCombiner, case05)
{
  //Test negative case: two not equal line numbers, not equal filenames, and equal tags
  std::string text1 =
    "This message does not matter because we do not filter by text";
  axom::lumberjack::Message m1;
  m1.text(text1);
  m1.addRank(13, 5);
  m1.fileName("foo.cpp");
  m1.lineNumber(12000);
  m1.tag("myTag");

  std::string text2 =
    "This message ALSO does not matter because we do not filter by text";
  axom::lumberjack::Message m2;
  m2.text(text2);
  m2.addRank(14, 5);
  m2.fileName("bar.cpp");
  m2.lineNumber(154);
  m2.tag("myTag");

  axom::lumberjack::LineFileTagCombiner c;

  bool shouldMessagesBeCombined = c.shouldMessagesBeCombined(m1, m2);

  if(shouldMessagesBeCombined == true)
  {
    c.combine(m1, m2, 5);
  }

  EXPECT_EQ(shouldMessagesBeCombined, false);

  EXPECT_EQ(m1.lineNumber(), 12000);
  EXPECT_EQ(m1.fileName().compare("foo.cpp"), 0);
  EXPECT_EQ(m1.text().compare(text1), 0);
  EXPECT_EQ(m1.count(), 1);
  EXPECT_EQ(m1.ranks().size(), 1);
  EXPECT_EQ(m1.ranks()[0], 13);
  EXPECT_EQ(m1.tag().compare("myTag"), 0);

  EXPECT_EQ(m2.lineNumber(), 154);
  EXPECT_EQ(m2.fileName().compare("bar.cpp"), 0);
  EXPECT_EQ(m2.text().compare(text2), 0);
  EXPECT_EQ(m2.count(), 1);
  EXPECT_EQ(m2.ranks().size(), 1);
  EXPECT_EQ(m2.ranks()[0], 14);
  EXPECT_EQ(m2.tag().compare("myTag"), 0);
}

TEST(lumberjack_LineFileTagCombiner, case06)
{
  //Test negative case: two not equal line numbers, equal filenames, and not equal tags
  std::string text1 =
    "This message does not matter because we do not filter by text";
  axom::lumberjack::Message m1;
  m1.text(text1);
  m1.addRank(13, 5);
  m1.fileName("foo.cpp");
  m1.lineNumber(12000);
  m1.tag("myTag");

  std::string text2 =
    "This message ALSO does not matter because we do not filter by text";
  axom::lumberjack::Message m2;
  m2.text(text2);
  m2.addRank(14, 5);
  m2.fileName("foo.cpp");
  m2.lineNumber(154);
  m2.tag("myOtherTag");

  axom::lumberjack::LineFileTagCombiner c;

  bool shouldMessagesBeCombined = c.shouldMessagesBeCombined(m1, m2);

  if(shouldMessagesBeCombined == true)
  {
    c.combine(m1, m2, 5);
  }

  EXPECT_EQ(shouldMessagesBeCombined, false);

  EXPECT_EQ(m1.lineNumber(), 12000);
  EXPECT_EQ(m1.fileName().compare("foo.cpp"), 0);
  EXPECT_EQ(m1.text().compare(text1), 0);
  EXPECT_EQ(m1.count(), 1);
  EXPECT_EQ(m1.ranks().size(), 1);
  EXPECT_EQ(m1.ranks()[0], 13);
  EXPECT_EQ(m1.tag().compare("myTag"), 0);

  EXPECT_EQ(m2.lineNumber(), 154);
  EXPECT_EQ(m2.fileName().compare("foo.cpp"), 0);
  EXPECT_EQ(m2.text().compare(text2), 0);
  EXPECT_EQ(m2.count(), 1);
  EXPECT_EQ(m2.ranks().size(), 1);
  EXPECT_EQ(m2.ranks()[0], 14);
  EXPECT_EQ(m2.tag().compare("myOtherTag"), 0);
}

TEST(lumberjack_LineFileTagCombiner, case07)
{
  //Test negative case: two equal line numbers, not equal filenames, and not equal tags
  std::string text1 =
    "This message does not matter because we do not filter by text";
  axom::lumberjack::Message m1;
  m1.text(text1);
  m1.addRank(13, 5);
  m1.fileName("foo.cpp");
  m1.lineNumber(154);
  m1.tag("myTag");

  std::string text2 =
    "This message ALSO does not matter because we do not filter by text";
  axom::lumberjack::Message m2;
  m2.text(text2);
  m2.addRank(14, 5);
  m2.fileName("bar.cpp");
  m2.lineNumber(154);
  m2.tag("myOtherTag");

  axom::lumberjack::LineFileTagCombiner c;

  bool shouldMessagesBeCombined = c.shouldMessagesBeCombined(m1, m2);

  if(shouldMessagesBeCombined == true)
  {
    c.combine(m1, m2, 5);
  }

  EXPECT_EQ(shouldMessagesBeCombined, false);

  EXPECT_EQ(m1.lineNumber(), 154);
  EXPECT_EQ(m1.fileName().compare("foo.cpp"), 0);
  EXPECT_EQ(m1.text().compare(text1), 0);
  EXPECT_EQ(m1.count(), 1);
  EXPECT_EQ(m1.ranks().size(), 1);
  EXPECT_EQ(m1.ranks()[0], 13);
  EXPECT_EQ(m1.tag().compare("myTag"), 0);

  EXPECT_EQ(m2.lineNumber(), 154);
  EXPECT_EQ(m2.fileName().compare("bar.cpp"), 0);
  EXPECT_EQ(m2.text().compare(text2), 0);
  EXPECT_EQ(m2.count(), 1);
  EXPECT_EQ(m2.ranks().size(), 1);
  EXPECT_EQ(m2.ranks()[0], 14);
  EXPECT_EQ(m2.tag().compare("myOtherTag"), 0);
}

TEST(lumberjack_LineFileTagCombiner, case08)
{
  //Test negative case: two not equal line numbers, not equal filenames, and not equal tags
  std::string text1 =
    "This message does not matter because we do not filter by text";
  axom::lumberjack::Message m1;
  m1.text(text1);
  m1.addRank(13, 5);
  m1.fileName("foo.cpp");
  m1.lineNumber(12000);
  m1.tag("myTag");

  std::string text2 =
    "This message ALSO does not matter because we do not filter by text";
  axom::lumberjack::Message m2;
  m2.text(text2);
  m2.addRank(14, 5);
  m2.fileName("bar.cpp");
  m2.lineNumber(154);
  m2.tag("myOtherTag");

  axom::lumberjack::LineFileTagCombiner c;

  bool shouldMessagesBeCombined = c.shouldMessagesBeCombined(m1, m2);

  if(shouldMessagesBeCombined == true)
  {
    c.combine(m1, m2, 5);
  }

  EXPECT_EQ(shouldMessagesBeCombined, false);

  EXPECT_EQ(m1.lineNumber(), 12000);
  EXPECT_EQ(m1.fileName().compare("foo.cpp"), 0);
  EXPECT_EQ(m1.text().compare(text1), 0);
  EXPECT_EQ(m1.count(), 1);
  EXPECT_EQ(m1.ranks().size(), 1);
  EXPECT_EQ(m1.ranks()[0], 13);
  EXPECT_EQ(m1.tag().compare("myTag"), 0);

  EXPECT_EQ(m2.lineNumber(), 154);
  EXPECT_EQ(m2.fileName().compare("bar.cpp"), 0);
  EXPECT_EQ(m2.text().compare(text2), 0);
  EXPECT_EQ(m2.count(), 1);
  EXPECT_EQ(m2.ranks().size(), 1);
  EXPECT_EQ(m2.ranks()[0], 14);
  EXPECT_EQ(m2.tag().compare("myOtherTag"), 0);
}