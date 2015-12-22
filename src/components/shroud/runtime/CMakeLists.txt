#
# Runtime library to support shroud.
#

#
# make the library
#
set(BUILD_SHARED_LIBS TRUE)
make_library(
    LIBRARY_NAME
        shroudrt
    LIBRARY_SOURCES
        shroudc.cpp
        fstr_mod.f
)

install(TARGETS shroudrt DESTINATION lib EXPORT ${PROJECT_NAME}-targets)
