#
# run-tutorial - A test program to use the tutorial examples
#

add_executable( shroud-run-tutorial
    main.f
    tutorial.cpp
)
target_link_libraries(shroud-run-tutorial shroudrt)
