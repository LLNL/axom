#
# Runtime library to support shroud.
#

#
# make the library
#
make_library( LIBRARY_NAME
                  shroudrt
              LIBRARY_SOURCES
                  fstr_mod.f
             )

install(TARGETS shroudrt DESTINATION lib EXPORT ${PROJECT_NAME}-targets)
