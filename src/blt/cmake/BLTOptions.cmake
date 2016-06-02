#
#
#

option(ENABLE_SHARED_LIBS "Enables shared libraries." OFF)

option(ENABLE_DOCS       "Enables documentation" ON)
option(ENABLE_EXAMPLES   "Enables examples" ON)
option(ENABLE_TESTS      "Enables tests" ON)
option(ENABLE_BENCHMARKS "Enables benchmarks" OFF)

option(ENABLE_COVERAGE "Enables code coverage support" OFF)
option(ENABLE_CXX11    "Enables C++11 language support" ON)
option(ENABLE_FORTRAN  "Enables Fortran compiler support" ON)

option(ENABLE_ALL_WARNINGS       "Enables all compiler warnings on all build targets" ON)
option(ENABLE_WARNINGS_AS_ERRORS "Enables treating compiler warnings as errors on all build targets" OFF)

option(ENABLE_MPI        "Enables MPI support" OFF)
option(ENABLE_OPENMP     "Enables OpenMP compiler support" OFF)

