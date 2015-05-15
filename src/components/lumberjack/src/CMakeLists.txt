configure_file(config.hpp.in "${HEADER_INCLUDES_DIRECTORY}/lumberjack/config.hpp")

#
# Specify all headers
#
set(lumberjack_headers
    Message.hpp
    )

#
# Specify all sources
#
set(lumberjack_sources
    Message.cpp
    )


#
# make the library
#
make_library(lumberjack "${lumberjack_sources}" "${lumberjack_headers}")
add_dependencies(lumberjack copy_headers_lumberjack)
install(TARGETS lumberjack DESTINATION lib EXPORT ${PROJECT_NAME}-targets)

#
# custom target to copy header files in the "include" folder of the build tree.
#
copy_headers_target(${PROJECT_NAME}  "${lumberjack_headers}"
                     ${HEADER_INCLUDES_DIRECTORY}/${PROJECT_NAME})
