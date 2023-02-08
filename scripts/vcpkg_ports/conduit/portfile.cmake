vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO llnl/conduit
    REF v0.8.6
    SHA512 b85c15bfa2687ba47f53c1ca269af72a1a31161848047e653bdc722a07f2682623640758cb5e83565ee655eca7cc993921c656208e6084513843927d76c5db66
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


