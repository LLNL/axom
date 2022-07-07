vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO mfem/mfem
    REF v4.4
    SHA512 aeb876de92638ae21b6d30c6e68d96bf0ad133a06664efd31a534e49915c71ce7cd214679f0e8021b557affac338c00444cca1b18fb35182a571a4967ae51009
    HEAD_REF master
    PATCHES "./export-extern-vars.patch"
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
else()
    # Move dll files from lib to bin directory
    file(MAKE_DIRECTORY ${CURRENT_PACKAGES_DIR}/bin )
    file(MAKE_DIRECTORY ${CURRENT_PACKAGES_DIR}/debug/bin )

    file(RENAME ${CURRENT_PACKAGES_DIR}/lib/mfem.dll
                ${CURRENT_PACKAGES_DIR}/bin/mfem.dll)

    file(RENAME ${CURRENT_PACKAGES_DIR}/debug/lib/mfem.dll
                ${CURRENT_PACKAGES_DIR}/debug/bin/mfem.dll)

    # Update paths to dlls in CMake config files
    foreach(_c  debug release)
        set(_f ${_config_dir}/MFEMTargets-${_c}.cmake)
        file(READ ${_f} _fdata)
        string(REPLACE "lib/mfem.dll" "bin/mfem.dll" _fdata "${_fdata}")
        file(WRITE  ${_f} "${_fdata}")
    endforeach()
endif()


# Put the license file where vcpkg expects it
file(INSTALL     ${SOURCE_PATH}/LICENSE 
     DESTINATION ${CURRENT_PACKAGES_DIR}/share/mfem 
     RENAME      copyright)


