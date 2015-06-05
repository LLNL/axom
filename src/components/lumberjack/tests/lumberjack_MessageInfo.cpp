#include "gtest/gtest.h"

#include "lumberjack/MessageInfo.hpp"

TEST(lumberjack_MessageInfo, getSet01)
{
	//Test most basic case: one message, one rank, file name, line number
	asctoolkit::lumberjack::MessageInfo m;
	m.message("I never wanted to do this job in the first place!");
	m.addRank(14, 5);
	m.fileName("foo.cpp");
	m.lineNumber(154);

	EXPECT_EQ(m.message(), "I never wanted to do this job in the first place!");
	EXPECT_EQ(m.fileName(), "foo.cpp");
	EXPECT_EQ(m.lineNumber(), 154);
	EXPECT_EQ(m.ranks().size(), (std::vector<int>::size_type)1);
	EXPECT_EQ(m.rankCount(), (std::vector<int>::size_type)1);
	EXPECT_EQ(m.ranks()[0], 14);
}

TEST(lumberjack_MessageInfo, getSet02)
{
	//Test that const char* will convert fine to string
	const char* messageConstCharPointer = "I... I wanted to be... A LUMBERJACK!";
	std::string messageString = "I... I wanted to be... A LUMBERJACK!";
	asctoolkit::lumberjack::MessageInfo m;
	m.message(messageConstCharPointer);
	m.addRank(14, 5);

	EXPECT_EQ(m.message(), messageString);
	EXPECT_EQ(m.ranks().size(), (std::vector<int>::size_type)1);
	EXPECT_EQ(m.rankCount(), (std::vector<int>::size_type)1);
	EXPECT_EQ(m.ranks()[0], 14);
}

TEST(lumberjack_MessageInfo, getSet03)
{
	//Test case: one message, filled ranks to rank limit
	const std::vector<int>::size_type rankLimit = 5;
	asctoolkit::lumberjack::MessageInfo m;
	m.message("Leaping from tree to tree! As they float down the mighty rivers of British Columbia!");
	for(int i=0; i<(int)rankLimit; ++i){
		m.addRank(i+1, rankLimit);
	}

	EXPECT_EQ(m.message(), "Leaping from tree to tree! As they float down the mighty rivers of British Columbia!");
	EXPECT_EQ(m.ranks().size(), rankLimit);
	EXPECT_EQ(m.rankCount(), rankLimit);
	for(int i=0; i<(int)rankLimit; ++i){
		EXPECT_EQ(m.ranks()[i], i+1);
	}
}

TEST(lumberjack_MessageInfo, getSet04)
{
	//Test case: one message, filled ranks to past rank limit
	const std::vector<int>::size_type rankLimit = 5;
	asctoolkit::lumberjack::MessageInfo m;
	m.message("With my best girl by my side!");
	for(int i=0; i<(int)rankLimit*2; ++i){
		m.addRank(i+1, rankLimit);
	}

	EXPECT_EQ(m.message(), "With my best girl by my side!");
	EXPECT_EQ(m.ranks().size(), rankLimit);
	EXPECT_EQ(m.rankCount(), rankLimit*2);
	for(int i=0; i<(int)rankLimit; ++i){
		EXPECT_EQ(m.ranks()[i], i+1);
	}
}

TEST(lumberjack_MessageInfo, getSet05)
{
	//Test case: one message, fill ranks with vector of 1 rank
	const std::vector<int>::size_type rankLimit = 5;
	std::vector<int> ranks;
	asctoolkit::lumberjack::MessageInfo m;
	m.message("The Larch! The Pine! The Giant Redwood tree! The Sequoia!");
	ranks.push_back(123);
	m.addRanks(ranks,rankLimit);

	EXPECT_EQ(m.message(), "The Larch! The Pine! The Giant Redwood tree! The Sequoia!");
	EXPECT_EQ(m.ranks().size(), (std::vector<int>::size_type)1);
	EXPECT_EQ(m.rankCount(), (std::vector<int>::size_type)1);
	EXPECT_EQ(m.ranks()[0], 123);
}

TEST(lumberjack_MessageInfo, getSet06)
{
	//Test case: one message, fill ranks with vector of ranks don't go past ranklimit
	const std::vector<int>::size_type rankLimit = 5;
	std::vector<int> ranks;
	asctoolkit::lumberjack::MessageInfo m;
	m.message("Oh, I'm a lumberjack, and I'm okay,");
	for(int i=0; i<(int)rankLimit; ++i){
		ranks.push_back(i+1);
	}
	m.addRanks(ranks,rankLimit);

	EXPECT_EQ(m.message(), "Oh, I'm a lumberjack, and I'm okay,");
	EXPECT_EQ(m.ranks().size(), rankLimit);
	EXPECT_EQ(m.rankCount(), rankLimit);
	for(int i=0; i<(int)rankLimit; ++i){
		EXPECT_EQ(m.ranks()[i], i+1);
	}
}

TEST(lumberjack_MessageInfo, getSet07)
{
	//Test case: one message, fill ranks with vector of ranks that will go past ranklimit
	const std::vector<int>::size_type rankLimit = 5;
	std::vector<int> ranks;
	asctoolkit::lumberjack::MessageInfo m;
	m.message("I sleep all night and I work all day.");
	for(int i=0; i<(int)rankLimit*3; ++i){
		ranks.push_back(i+1);
	}
	m.addRanks(ranks,rankLimit);

	EXPECT_EQ(m.message(), "I sleep all night and I work all day.");
	EXPECT_EQ(m.ranks().size(), rankLimit);
	EXPECT_EQ(m.rankCount(), rankLimit*3);
	for(int i=0; i<(int)rankLimit; ++i){
		EXPECT_EQ(m.ranks()[i], i+1);
	}
}

TEST(lumberjack_MessageInfo, testConstructor01)
{
	//Test most basic case: one message, one rank, file name, line number
	asctoolkit::lumberjack::MessageInfo m("He's a lumberjack, and he's okay,",
									  122, "foo.cpp", 154);

	EXPECT_EQ(m.message(), "He's a lumberjack, and he's okay,");
	EXPECT_EQ(m.fileName(), "foo.cpp");
	EXPECT_EQ(m.lineNumber(), 154);
	EXPECT_EQ(m.ranks().size(), (std::vector<int>::size_type)1);
	EXPECT_EQ(m.rankCount(), (std::vector<int>::size_type)1);
	EXPECT_EQ(m.ranks()[0], 122);
}

TEST(lumberjack_MessageInfo, testConstructor02)
{
	//Test most basic case: one message, one rank, file name, line number
	const std::vector<int>::size_type rankLimit = 5;
	std::vector<int> ranks;
	for(int i=0; i<(int)rankLimit; ++i){
		ranks.push_back(i+1);
	}

	asctoolkit::lumberjack::MessageInfo m("He sleeps all night and he works all day.",
									  ranks, "foo.cpp", 154);

	EXPECT_EQ(m.message(), "He sleeps all night and he works all day.");
	EXPECT_EQ(m.fileName(), "foo.cpp");
	EXPECT_EQ(m.lineNumber(), 154);
	EXPECT_EQ(m.ranks().size(), rankLimit);
	EXPECT_EQ(m.rankCount(), rankLimit);
	for(int i=0; i<(int)rankLimit; ++i){
		EXPECT_EQ(m.ranks()[i], i+1);
	}
}
