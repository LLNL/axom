################
# C++ Benchmarks
################
set(slam_benchmark_files
    slam_array.cpp
    slam_sets.cpp
    )

foreach(test ${slam_benchmark_files})
  blt_add_benchmark( TEST_SOURCE ${test} DEPENDS_ON slam slic boost)
endforeach()

