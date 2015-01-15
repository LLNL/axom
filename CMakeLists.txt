cmake_minimum_required(VERSION 2.8)

#SET(CMAKE_CXX_COMPILER  /usr/apps/gnu/4.9.0/bin/g++)
SET(CMAKE_CXX_COMPILER  /usr/local/bin/g++4.9.2)

project(test)

add_executable( ds_test src/DataGroup.cpp  src/DataObject.cpp src/DatastoreInterface.cpp src/Examples/ShockTube.cc src/test.cpp)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g3 -std=c++11 -Wall -Wabi -Wextra -Wshadow -pedantic  -Woverloaded-virtual  -Wno-unused-variable")
SET(CMAKE_VERBOSE_MAKEFILE)
include_directories(src)

install(TARGETS ds_test DESTINATION bin)
