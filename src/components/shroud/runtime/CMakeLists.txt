#
# Runtime library to support shroud.
#
set(shroud_headers
    shroudrt.hpp
)

#
# custom target to copy header files in the "include" folder of the build tree.
#
#copy_headers_target( ${PROJECT_NAME} "${shroud_headers}"
#                     ${HEADER_INCLUDES_DIRECTORY}/${PROJECT_NAME} )
#set(PROJECT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
#copy_headers_target( shroudrt "${shroud_headers}"
#                     ${HEADER_INCLUDES_DIRECTORY}/shroud )

#
# make the library
#
make_library(
    LIBRARY_NAME
        shroudrt
    LIBRARY_SOURCES
        ${shroud_headers}
        shroudc.cpp
        fstr_mod.f
)

install(TARGETS shroudrt DESTINATION lib EXPORT ${PROJECT_NAME}-targets)
