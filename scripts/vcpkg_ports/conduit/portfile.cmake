if(VCPKG_CMAKE_SYSTEM_NAME STREQUAL "WindowsStore")
    message(FATAL_ERROR "${PORT} does not currently support UWP")
endif()

include(vcpkg_common_functions)
vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO llnl/conduit
    REF v0.6.0
    SHA512 b65ce01a5aabee660b62bf80daf629af200e6c489bdb2e5c359154bc212c740cca663d421c12181b2d8fe40a694430f246f3fc052c6945b430892e01e58b59f0
    HEAD_REF master
)

set(_is_shared TRUE)
if(VCPKG_LIBRARY_LINKAGE STREQUAL "static")
    set(_is_shared FALSE)
endif()

vcpkg_configure_cmake(
    SOURCE_PATH ${SOURCE_PATH}/src
    PREFER_NINJA
    OPTIONS 
        -DBLT_SOURCE_DIR=${CURRENT_INSTALLED_DIR}/share/blt
        -DCONDUIT_ENABLE_TESTS=OFF
        -DENABLE_COVERAGE=OFF
        -DENABLE_DOCS=OFF
        -DENABLE_EXAMPLES=OFF
        -DENABLE_PYTHON=OFF
        -DENABLE_TESTS=OFF
        -DENABLE_UTILS=OFF
        -DBUILD_SHARED_LIBS=${_is_shared}
        -DHDF5_DIR=${CURRENT_INSTALLED_DIR}
)

vcpkg_install_cmake()
vcpkg_fixup_cmake_targets(CONFIG_PATH lib/cmake)
vcpkg_copy_pdbs()


## shuffle the output directories to make vcpkg happy
# Remove extraneous debug header files
file(REMOVE_RECURSE ${CURRENT_PACKAGES_DIR}/debug/include)
file(REMOVE_RECURSE ${CURRENT_PACKAGES_DIR}/debug/share)

# Remove exe files -- vcpkg doesn't like them
# (Future): It might be possible to move them to the vcpkg 'tools' directory
foreach(_dir "bin" "debug/bin")
    file(GLOB _exe ${CURRENT_PACKAGES_DIR}/${_dir}/*.exe)
    if(_exe)
        file(REMOVE ${_exe})
    endif()
endforeach()

if(VCPKG_LIBRARY_LINKAGE STREQUAL static)
    file(REMOVE_RECURSE ${CURRENT_PACKAGES_DIR}/bin ${CURRENT_PACKAGES_DIR}/debug/bin)
endif()


# Move/shuffle cmake files to 'share' directory
#file(MAKE_DIRECTORY ${CURRENT_PACKAGES_DIR}/share/conduit )
#file(COPY ${CURRENT_PACKAGES_DIR}/lib/cmake/ 
#     DESTINATION ${CURRENT_PACKAGES_DIR}/share/conduit/ ) 
#file(RENAME ${CURRENT_PACKAGES_DIR}/debug/lib/cmake/conduit-debug.cmake 
#            ${CURRENT_PACKAGES_DIR}/share/conduit/)
#file(REMOVE_RECURSE ${CURRENT_PACKAGES_DIR}/lib/cmake)
#file(REMOVE_RECURSE ${CURRENT_PACKAGES_DIR}/debug/lib/cmake)

# TODO: Fixup cmake files or config files?

# Put the license file where vcpkg expects it
file(INSTALL ${SOURCE_PATH}/LICENSE DESTINATION ${CURRENT_PACKAGES_DIR}/share/conduit RENAME copyright)


