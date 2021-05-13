// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "benchmark/benchmark_api.h"
#include "axom/slic/interface/slic.hpp"
#include "axom/slic/core/SimpleLogger.hpp"
using axom::slic::SimpleLogger;

/*!
 * \file
 *
 * A series of tests of SLIC logging macros (SLIC_WARNING to avoid killing the
 * program) in the TPL benchmark framework.
 *
 * These tests were written to see why SLIC macros was not being logged
 * when called from a benchmark fixture's destructor.
 * It turns out that the fixtures are destructed after main exits,
 * so SLIC is already finalized.
 *
 */

namespace
{
static const bool IS_VALID = true;
static const bool SHOULD_PRINT = false;

void printMsg(std::string const& str, bool override = false)
{
  if(SHOULD_PRINT || override) std::cout << str << std::endl;
}

struct AssertCtor
{
  AssertCtor() { SLIC_WARNING_IF(IS_VALID, "Testing warning in .ctor"); }
};

struct AssertMethod
{
  void foo() { SLIC_WARNING_IF(IS_VALID, "Testing warning in class method"); }
};

struct AssertDtor
{
  ~AssertDtor() { SLIC_WARNING_IF(IS_VALID, "Testing warning in .dtor"); }
};

class SetFixtureC : public ::benchmark::Fixture
{
public:
  SetFixtureC()
  {
    SLIC_WARNING_IF(IS_VALID, "Testing warning in .ctor");
    printMsg(
      "*** Note: the slic log is not printed for SetFixtureC/assertCtor since "
      "it is called before slic::initialize() ***",
      true);
  }
};

class SetFixtureS : public ::benchmark::Fixture
{
public:
  void SetUp() { SLIC_WARNING_IF(IS_VALID, "Testing warning in setup"); }
};

class SetFixtureT : public ::benchmark::Fixture
{
public:
  void TearDown() { SLIC_WARNING_IF(IS_VALID, "Testing warning in teardown"); }
};

class SetFixtureD : public ::benchmark::Fixture
{
public:
  ~SetFixtureD()
  {
    SLIC_WARNING_IF(IS_VALID, "Testing warning in .dtor");
    printMsg(
      "*** Note: the slic log is not printed for SetFixtureD/assertDtor since "
      "it is called after slic::finalize() ***",
      true);
  }
};

class SetFixtureOutput : public ::benchmark::Fixture
{
public:
  SetFixtureOutput() { printMsg("  fixure .ctor"); }
  void SetUp() { printMsg("  fixure SetUp"); }
  void TearDown() { printMsg("  fixure TearDown"); }
  ~SetFixtureOutput() { printMsg("  fixure ~.dtor"); }
};

static void BM_AssertBefore(benchmark::State& state)
{
  SLIC_WARNING_IF(IS_VALID, "Testing warning before running benchmark");
  while(state.KeepRunning()) benchmark::DoNotOptimize(0);
}

static void BM_AssertAfter(benchmark::State& state)
{
  while(state.KeepRunning()) benchmark::DoNotOptimize(0);
  SLIC_WARNING_IF(IS_VALID, "Testing warning after running benchmark");
}

static void BM_CallAssertCtorBefore(benchmark::State& state)
{
  AssertCtor();
  while(state.KeepRunning()) benchmark::DoNotOptimize(0);
}

static void BM_CallAssertDtorBefore(benchmark::State& state)
{
  AssertDtor();
  while(state.KeepRunning()) benchmark::DoNotOptimize(0);
}

static void BM_CallAssertDtorDuring(benchmark::State& state)
{
  while(state.KeepRunning()) AssertDtor();
}

}  // namespace

BENCHMARK(BM_AssertBefore);
BENCHMARK(BM_AssertAfter);
BENCHMARK(BM_CallAssertCtorBefore);
BENCHMARK(BM_CallAssertDtorBefore);
BENCHMARK(BM_CallAssertDtorDuring);

BENCHMARK_F(SetFixtureC, assertCtor)(benchmark::State& state)
{
  while(state.KeepRunning())
  {
    benchmark::DoNotOptimize(0);
  }
}
BENCHMARK_F(SetFixtureS, assertSetup)(benchmark::State& state)
{
  while(state.KeepRunning())
  {
    benchmark::DoNotOptimize(0);
  }
}
BENCHMARK_F(SetFixtureT, assertTeardown)(benchmark::State& state)
{
  while(state.KeepRunning())
  {
    benchmark::DoNotOptimize(0);
  }
}
BENCHMARK_F(SetFixtureD, assertDtor)(benchmark::State& state)
{
  while(state.KeepRunning())
  {
    benchmark::DoNotOptimize(0);
  }
}

// The following two benchmarks are here to show the order of construction,
// setup, teardown and destruction of benchmarks fixtures.
// Note: Need to enable printing (by setting SHOULD_PRINT above to true) to see
// this.
BENCHMARK_F(SetFixtureOutput, process1)(benchmark::State& state)
{
  printMsg("   P1 (before)");
  while(state.KeepRunning())
  {
    benchmark::DoNotOptimize(0);
  }
  printMsg("   P1 (after)");
}
BENCHMARK_F(SetFixtureOutput, process2)(benchmark::State& state)
{
  printMsg("   P2 (before)");
  while(state.KeepRunning())
  {
    benchmark::DoNotOptimize(0);
  }
  printMsg("   P2 (after)");
}

int main(int argc, char* argv[])
{
  printMsg("at start of main");

  printMsg(" Before init logger");
  axom::slic::SimpleLogger logger;  // create & initialize test logger,
  printMsg(" After init logger");

  printMsg(" Before init benchmark");

  ::benchmark::Initialize(&argc, argv);
  printMsg(" after init benchmark");

  printMsg(" before run benchmark");
  ::benchmark::RunSpecifiedBenchmarks();
  printMsg(" after run benchmark");

  printMsg("about to return from main");
  return 0;
}
