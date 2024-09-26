vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO mfem/mfem
    REF v4.6
    SHA512 8805b4993b6f11abe7ac7dda59d0ddb2e0f5f6b09c2b9c57e665f481cd9bd6b669e63621b38989f70dc8ae38c42a7e8c4e10a1d87a4ac29d53ddd95ce79db0ae
    HEAD_REF master
    )

set(_is_shared TRUE)
if(VCPKG_LIBRARY_LINKAGE STREQUAL "static")
    set(_is_shared FALSE)
endif()

vcpkg_configure_cmake(
    SOURCE_PATH ${SOURCE_PATH}
    PREFER_NINJA
    OPTIONS 
        -DMFEM_ENABLE_EXAMPLES=OFF
        -DMFEM_ENABLE_MINIAPPS=OFF
        -DMFEM_ENABLE_TESTING=OFF
        -DBUILD_SHARED_LIBS=${_is_shared}
        -DCMAKE_WINDOWS_EXPORT_ALL_SYMBOLS=${_is_shared}
)

vcpkg_install_cmake()
vcpkg_fixup_cmake_targets(CONFIG_PATH lib/cmake)
vcpkg_copy_pdbs()


## shuffle the output directories to make vcpkg happy
# Remove extraneous debug header files
file(REMOVE_RECURSE ${CURRENT_PACKAGES_DIR}/debug/include)
file(REMOVE_RECURSE ${CURRENT_PACKAGES_DIR}/debug/share)

# Move CMake config files up a directory
set(_config_dir "${CURRENT_PACKAGES_DIR}/share/mfem")
file(GLOB _cmake_files "${_config_dir}/mfem/*.cmake")
foreach(_f ${_cmake_files})
    get_filename_component(_name ${_f} NAME)
    file(RENAME ${_f} ${_config_dir}/${_name})
endforeach()
file(REMOVE_RECURSE "${_config_dir}/mfem")

if(VCPKG_LIBRARY_LINKAGE STREQUAL static)
    # Note: Not tested
    file(REMOVE_RECURSE ${CURRENT_PACKAGES_DIR}/bin ${CURRENT_PACKAGES_DIR}/debug/bin)
endif()


# Put the license file where vcpkg expects it
file(INSTALL     ${SOURCE_PATH}/LICENSE 
     DESTINATION ${CURRENT_PACKAGES_DIR}/share/mfem 
     RENAME      copyright)


