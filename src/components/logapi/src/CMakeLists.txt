#
# Specify all headers
#
set(logapi_headers
    Console.hpp
    Logger.hpp
    LogStream.hpp
    )

#
# Specify all sources
#
set(logapi_sources
    Console.cpp
    Logger.cpp
    LogStream.cpp
    )

if (ENABLE_MPI)
   list(APPEND logapi_headers SynchronizedConsole.hpp)
   list(APPEND logapi_sources SynchronizedConsole.cpp)
endif()

#
# make the library
#
make_library(logapi "${logapi_sources}" "${logapi_headers}")
add_dependencies(logapi copy_headers_logapi)
install( TARGETS logapi DESTINATION lib EXPORT ${PROJECT_NAME}-targets)

#
# custom target to copy header files in the "include" folder of the build tree.
#
copy_headers_target( ${PROJECT_NAME}
                     "${logapi_headers}"
                     ${HEADER_INCLUDES_DIRECTORY}/${PROJECT_NAME} )