if(VCPKG_CMAKE_SYSTEM_NAME STREQUAL "WindowsStore")
    message(FATAL_ERROR "${PORT} does not currently support UWP")
endif()

include(vcpkg_common_functions)
vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO llnl/conduit
    REF v0.4.0
    SHA512 772c932ef2265ba244bd06feb29f8d81bc64920e7df581818d9292a6a6007428b6f883bfb3810dba9402f0a38220eae634c75ab966b7c1d96fcd63f3f5f5589f
    HEAD_REF master
)

message(STATUS "SOURCE_PATH -- ${SOURCE_PATH}")

# TODO: Needs to clone recursively for BLT!

vcpkg_configure_cmake(
    SOURCE_PATH ${SOURCE_PATH}/src
)

vcpkg_install_cmake()
