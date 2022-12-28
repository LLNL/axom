vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO llnl/camp
    REF v2022.03.2
    SHA512 d124c0e933f042525e9b48c21b93e7f4f634dfc0f87742e018a1c7de43ed6e30670d8f8e4ce369018a8e1c884b2b27f4755ee9f07a077820b2a3133546f6d622
    HEAD_REF develop
)

set(_is_shared TRUE)
if(VCPKG_LIBRARY_LINKAGE STREQUAL "static")
    set(_is_shared FALSE)
endif()

vcpkg_configure_cmake(
    SOURCE_PATH ${SOURCE_PATH}
    PREFER_NINJA
    OPTIONS 
        -DBLT_SOURCE_DIR=${CURRENT_INSTALLED_DIR}/share/blt
        -DENABLE_COVERAGE=OFF
        -DENABLE_ALL_WARNINGS=OFF
        -DENABLE_DOCS=OFF
        -DENABLE_EXAMPLES=OFF
        -DENABLE_TESTS=OFF
        -DBUILD_SHARED_LIBS=${_is_shared}
)

vcpkg_install_cmake()
vcpkg_cmake_config_fixup()
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

# Put the license file where vcpkg expects it
file(INSTALL ${SOURCE_PATH}/LICENSE DESTINATION ${CURRENT_PACKAGES_DIR}/share/camp RENAME copyright)
