/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */


#include <cstdlib>
#include <ctime>

#include "benchmark/benchmark_api.h"
#include "slic/slic.hpp"
#include "slic/UnitTestLogger.hpp"

//------------------------------------------------------------------------------
namespace {
    const int STRIDE = 7;
    const int OFFSET = 12;


    typedef int DataType;
    typedef DataType* DataArray;

    DataArray generateRandomArray(int sz)
    {
        DataArray data = new DataType[sz];
        for(int i=0; i< sz; ++i)
        {
            data[i] = rand();
        }

        return data;
    }

}
//------------------------------------------------------------------------------


#define BASIC_BENCHMARK_TEST(x) \
    BENCHMARK(x)->Arg( 1<<3 )->Arg( 1<<9 )->Arg( 1 << 13 )->Arg( 1 << 24 )


void benchmark_contig_sequence(benchmark::State& state) {
    while (state.KeepRunning()) {
        for (int i=0; i < state.range_x(); ++i) {
          benchmark::DoNotOptimize(i);
      }
  }
  state.SetItemsProcessed(state.iterations() * state.range_x());

}
BASIC_BENCHMARK_TEST(benchmark_contig_sequence);

void benchmark_strided_sequence(benchmark::State& state) {
    while (state.KeepRunning()) {
        for (int i=0; i < state.range_x(); i += STRIDE) {
          benchmark::DoNotOptimize(i);
      }
  }
  state.SetItemsProcessed(state.iterations() * state.range_x());

}
BASIC_BENCHMARK_TEST(benchmark_strided_sequence);


void benchmark_offset_sequence(benchmark::State& state) {
    while (state.KeepRunning()) {
        for (int i=OFFSET; i < state.range_x(); ++i) {
          benchmark::DoNotOptimize(i);
      }
  }
  state.SetItemsProcessed(state.iterations() * state.range_x());

}
BASIC_BENCHMARK_TEST(benchmark_offset_sequence);


void benchmark_offset_strided_sequence(benchmark::State& state) {
    while (state.KeepRunning()) {
        for (int i=OFFSET; i < state.range_x(); i+= STRIDE) {
          benchmark::DoNotOptimize(i);
      }
  }
  state.SetItemsProcessed(state.iterations() * state.range_x());

}
BASIC_BENCHMARK_TEST(benchmark_offset_strided_sequence);


void benchmark_indirection_sequence(benchmark::State& state)
{
    DataArray data = generateRandomArray( state.range_x() );

    while (state.KeepRunning())
    {

        for (int i=0; i < state.range_x(); ++i)
        {
          benchmark::DoNotOptimize( data[i] );
        }

    }

    delete[] data;

  state.SetItemsProcessed(state.iterations() * state.range_x());

}
BASIC_BENCHMARK_TEST(benchmark_indirection_sequence);







int main(int argc, char * argv[])
{
  std::srand (std::time(NULL));
  asctoolkit::slic::UnitTestLogger logger;  // create & initialize test logger,

  ::benchmark::Initialize(&argc, argv);
  ::benchmark::RunSpecifiedBenchmarks();

  return 0;
}
