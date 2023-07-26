vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO llnl/umpire
    REF v2022.10.0
    SHA512 38bc74c360ee73e8ee04fbb652cff160551597a46af587f8c9da755cef591ae4add66d9af038a0a14722653d135007b01b37d3addf4d64ca0d1ed129e0461428
    HEAD_REF develop
)

vcpkg_check_features(OUT_FEATURE_OPTIONS FEATURE_OPTIONS
    FEATURES
        openmp       ENABLE_OPENMP
)

set(_is_shared TRUE)
if(VCPKG_LIBRARY_LINKAGE STREQUAL "static")
    set(_is_shared FALSE)
else()
    list(APPEND FEATURE_OPTIONS -DCMAKE_WINDOWS_EXPORT_ALL_SYMBOLS=ON)
endif()


vcpkg_configure_cmake(
    SOURCE_PATH ${SOURCE_PATH}
    PREFER_NINJA
    OPTIONS 
        -DBLT_SOURCE_DIR:PATH=${CURRENT_INSTALLED_DIR}/share/blt
        -Dcamp_DIR:PATH=${CURRENT_INSTALLED_DIR}
        -DENABLE_ALL_WARNINGS:BOOL=OFF
        -DENABLE_WARNINGS_AS_ERRORS:BOOL=OFF
        -DENABLE_COVERAGE:BOOL=OFF
        -DENABLE_EXAMPLES:BOOL=OFF
        -DENABLE_TESTS:BOOL=OFF
        -DENABLE_BENCHMARKS:BOOL=OFF
        -DUMPIRE_ENABLE_FILESYSTEM:BOOL=ON
        -DBLT_CXX_STD:STRING=c++17
        -DBLT_OPENMP_LINK_FLAGS:STRING=" "
        -DUMPIRE_ENABLE_TOOLS:BOOL=OFF
        -DUMPIRE_ENABLE_TESTS:BOOL=OFF
        -DUMPIRE_ENABLE_BENCHMARKS:BOOL=OFF
        -DUMPIRE_ENABLE_DEVELOPER_BENCHMARKS:BOOL=OFF
        -DBUILD_SHARED_LIBS:BOOL=${_is_shared}
        ${FEATURE_OPTIONS}
)

vcpkg_install_cmake()
vcpkg_fixup_cmake_targets(CONFIG_PATH lib/cmake/umpire
                          TARGET_PATH share/umpire)
vcpkg_copy_pdbs()


## shuffle the output directories to make vcpkg happy
# Remove extraneous debug header files
file(REMOVE_RECURSE ${CURRENT_PACKAGES_DIR}/debug/include)
file(REMOVE_RECURSE ${CURRENT_PACKAGES_DIR}/debug/share)

message(STATUS "CURRENT_PACKAGES_DIR -- ${CURRENT_PACKAGES_DIR}")

if(VCPKG_LIBRARY_LINKAGE STREQUAL static)
    # Note: Not tested
    file(REMOVE_RECURSE ${CURRENT_PACKAGES_DIR}/bin ${CURRENT_PACKAGES_DIR}/debug/bin)
else()
    set(_config_dir "${CURRENT_PACKAGES_DIR}/share/umpire")
    
    # Move dll files from lib to bin directory
    file(MAKE_DIRECTORY ${CURRENT_PACKAGES_DIR}/bin )
    file(MAKE_DIRECTORY ${CURRENT_PACKAGES_DIR}/debug/bin )

    file(RENAME ${CURRENT_PACKAGES_DIR}/lib/umpire.dll
                ${CURRENT_PACKAGES_DIR}/bin/umpire.dll)

    file(RENAME ${CURRENT_PACKAGES_DIR}/debug/lib/umpire.dll
                ${CURRENT_PACKAGES_DIR}/debug/bin/umpire.dll)

    # Update paths to dlls in CMake config files
    foreach(_c  debug release)
        set(_f ${_config_dir}/umpire-targets-${_c}.cmake)
        file(READ ${_f} _fdata)
        string(REPLACE "lib/umpire.dll" "bin/umpire.dll" _fdata "${_fdata}")
        file(WRITE  ${_f} "${_fdata}")
    endforeach()

    # Fix erroneous "include" path appended as system include directories
    set(_f ${_config_dir}/umpire-targets.cmake)
    set(_pattern "INTERFACE_SYSTEM_INCLUDE_DIRECTORIES \"include.*")
    file(STRINGS ${_f} _file_lines)
    file(WRITE ${_f} "")
    foreach(_line IN LISTS _file_lines)
        #message(STATUS "\t\tRegex input:  ${_line}")
        string(REGEX REPLACE ${_pattern} "" _stripped "${_line}")
        #message(STATUS "\t\tRegex output: ${_stripped}")
        file(APPEND ${_f} "${_stripped}\n")
    endforeach()
endif()


# Put the license file where vcpkg expects it
file(INSTALL ${SOURCE_PATH}/LICENSE DESTINATION ${CURRENT_PACKAGES_DIR}/share/umpire RENAME copyright)
