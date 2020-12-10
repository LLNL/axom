// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/Units.hpp"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <stdexcept>

namespace axom
{
namespace klee
{
using ::testing::DoubleEq;
using ::testing::ElementsAre;
using ::testing::Matches;

struct Length
{
  double value;
  LengthUnit units;
};

bool convertsTo(Length length1, Length length2)
{
  double length2InLength1Units =
    convert(length2.value, length2.units, length1.units);
  return Matches(DoubleEq(length1.value))(length2InLength1Units);
}

bool areEquivalent(Length length1, Length length2)
{
  return convertsTo(length1, length2) && convertsTo(length2, length1);
}

TEST(Units, parseLengthUnits)
{
  EXPECT_EQ(LengthUnit::km, parseLengthUnits("km"));
  EXPECT_EQ(LengthUnit::m, parseLengthUnits("m"));
  EXPECT_EQ(LengthUnit::dm, parseLengthUnits("dm"));
  EXPECT_EQ(LengthUnit::cm, parseLengthUnits("cm"));
  EXPECT_EQ(LengthUnit::mm, parseLengthUnits("mm"));
  EXPECT_EQ(LengthUnit::um, parseLengthUnits("um"));
  EXPECT_EQ(LengthUnit::nm, parseLengthUnits("nm"));
  EXPECT_EQ(LengthUnit::angstrom, parseLengthUnits("A"));
  EXPECT_EQ(LengthUnit::miles, parseLengthUnits("miles"));
  EXPECT_EQ(LengthUnit::feet, parseLengthUnits("ft"));
  EXPECT_EQ(LengthUnit::feet, parseLengthUnits("feet"));
  EXPECT_EQ(LengthUnit::inches, parseLengthUnits("in"));
  EXPECT_EQ(LengthUnit::inches, parseLengthUnits("inches"));
  EXPECT_EQ(LengthUnit::mils, parseLengthUnits("mils"));

  EXPECT_THROW(parseLengthUnits("no such units"), std::invalid_argument);
}

TEST(Units, getConversionFactor)
{
  EXPECT_DOUBLE_EQ(1000, getConversionFactor(LengthUnit::km, LengthUnit::m));
  EXPECT_DOUBLE_EQ(10, getConversionFactor(LengthUnit::m, LengthUnit::dm));
  EXPECT_DOUBLE_EQ(10, getConversionFactor(LengthUnit::dm, LengthUnit::cm));
  EXPECT_DOUBLE_EQ(10, getConversionFactor(LengthUnit::cm, LengthUnit::mm));
  EXPECT_DOUBLE_EQ(1000, getConversionFactor(LengthUnit::mm, LengthUnit::um));
  EXPECT_DOUBLE_EQ(1000, getConversionFactor(LengthUnit::um, LengthUnit::nm));
  EXPECT_DOUBLE_EQ(10, getConversionFactor(LengthUnit::nm, LengthUnit::angstrom));
  EXPECT_DOUBLE_EQ(2.54, getConversionFactor(LengthUnit::inches, LengthUnit::cm));
  EXPECT_DOUBLE_EQ(1000,
                   getConversionFactor(LengthUnit::inches, LengthUnit::mils));
  EXPECT_DOUBLE_EQ(12, getConversionFactor(LengthUnit::feet, LengthUnit::inches));
  EXPECT_DOUBLE_EQ(5280,
                   getConversionFactor(LengthUnit::miles, LengthUnit::feet));

  EXPECT_DOUBLE_EQ(1, getConversionFactor(LengthUnit::miles, LengthUnit::miles));

  EXPECT_THROW(getConversionFactor(LengthUnit::cm, LengthUnit::unspecified),
          std::invalid_argument);
  EXPECT_THROW(getConversionFactor(LengthUnit::unspecified, LengthUnit::cm),
          std::invalid_argument);
}

TEST(Units, convert_adjacent)
{
  EXPECT_TRUE(areEquivalent({1, LengthUnit::km}, {1000, LengthUnit::m}));
  EXPECT_TRUE(areEquivalent({1, LengthUnit::m}, {10, LengthUnit::dm}));
  EXPECT_TRUE(areEquivalent({1, LengthUnit::dm}, {10, LengthUnit::cm}));
  EXPECT_TRUE(areEquivalent({1, LengthUnit::cm}, {10, LengthUnit::mm}));
  EXPECT_TRUE(areEquivalent({1, LengthUnit::mm}, {1000, LengthUnit::um}));
  EXPECT_TRUE(areEquivalent({1, LengthUnit::um}, {1000, LengthUnit::nm}));
  EXPECT_TRUE(areEquivalent({1, LengthUnit::nm}, {10, LengthUnit::angstrom}));
  EXPECT_TRUE(areEquivalent({1, LengthUnit::inches}, {2.54, LengthUnit::cm}));
  EXPECT_TRUE(areEquivalent({1, LengthUnit::inches}, {1000, LengthUnit::mils}));
  EXPECT_TRUE(areEquivalent({1, LengthUnit::feet}, {12, LengthUnit::inches}));
  EXPECT_TRUE(areEquivalent({1, LengthUnit::miles}, {5280, LengthUnit::feet}));
}

TEST(Units, convert_all)
{
  Length equivalentLengths[] = {
    {0.254254e-3, LengthUnit::km},
    {0.254254, LengthUnit::m},
    {0.254254e1, LengthUnit::dm},
    {0.254254e2, LengthUnit::cm},
    {0.254254e3, LengthUnit::mm},
    {0.254254e6, LengthUnit::um},
    {0.254254e9, LengthUnit::nm},
    {0.254254e10, LengthUnit::angstrom},
    {10.01, LengthUnit::inches},
    {10.01e3, LengthUnit::mils},
    {10.01 / 12, LengthUnit::feet},
    {10.01 / 12 / 5280, LengthUnit::miles},
  };

  for(auto source : equivalentLengths)
  {
    for(auto target : equivalentLengths)
    {
      EXPECT_TRUE(areEquivalent(source, target));
    }
  }
}

TEST(Units, convertAll)
{
  std::vector<double> values {1.0, 2.0, 3.0};
  convertAll(values, LengthUnit::m, LengthUnit::cm);
  EXPECT_THAT(values, ElementsAre(100, 200, 300));
}

}  // namespace klee
}  // namespace axom