configure_file(config.hpp.in "${HEADER_INCLUDES_DIRECTORY}/lumberjack/config.hpp")

#
# Specify all headers
#
set(lumberjack_headers
	Combiner.hpp
	Communicator.hpp
	MessageEqualityCombiner.hpp
    MessageInfo.hpp
    RootCommunicator.hpp
    )

#
# Specify all sources
#
set(lumberjack_sources
    MessageInfo.cpp
    RootCommunicator.cpp
    )


#
# make the library
#
make_library( LIBRARY_NAME
                  lumberjack
              LIBRARY_SOURCES
                  "${lumberjack_sources}"
                  "${lumberjack_headers}"
              WITH_MPI
              )

add_dependencies(lumberjack copy_headers_lumberjack)
add_dependencies(lumberjack copy_headers_common)

install(TARGETS lumberjack DESTINATION lib EXPORT ${PROJECT_NAME}-targets)

#
# custom target to copy header files in the "include" folder of the build tree.
#
copy_headers_target(${PROJECT_NAME}  "${lumberjack_headers}"
                     ${HEADER_INCLUDES_DIRECTORY}/${PROJECT_NAME})
