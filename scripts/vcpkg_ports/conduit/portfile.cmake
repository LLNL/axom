vcpkg_fail_port_install(ON_TARGET "UWP")

include(vcpkg_common_functions)
vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO llnl/conduit
    REF v0.7.2
    SHA512 85ee7c2a5483952c0f8b093a0fef95f393cd5c34884caefb03734cfd4095b8f64440b8c599a0547a25adea6da5b0126445885a681d706d4b5b68c22d897aea18
    HEAD_REF develop
    PATCHES "./hdf5-1-12.patch"
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


