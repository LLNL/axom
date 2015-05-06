## Build examples using logapi
make_executable( EXECUTABLE_SOURCE logging_example.cc DEPENDS_ON logapi )

## Build MPI examples using logapi
if (ENABLE_MPI)
   make_executable( EXECUTABLE_SOURCE
                        parallel_logging_example.cc
                    DEPENDS_ON
                        logapi ${MPI_C_LIBRARIES}  )
endif()