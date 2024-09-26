vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO llnl/conduit
    REF v0.9.2
    SHA512 f067423373ad542ba26c6331e17f78d10ad790f33b68c4c44dc54839ab1ef7e68bec5580de18ac07d5609a745f5b285822b6061e7c144e5a7d487864a605a874
    HEAD_REF develop
    PATCHES 
        "./setup_deps_vcpkg_triplet.patch"
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
        -DCONDUIT_INSTALL_CONFIG_DIR="share/conduit"
        -DCONDUIT_INSTALL_CMAKE_MODULE_DIR="share"
)

vcpkg_install_cmake()
vcpkg_fixup_cmake_targets(CONFIG_PATH lib/cmake
                          TARGET_PATH share
                          TOOLS_PATH tools/conduit)
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


