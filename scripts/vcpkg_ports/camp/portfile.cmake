vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO llnl/camp
    REF v2025.03.0
    SHA512 8d1b0f89a4ca95d418a6bf0b9df7e417b9130b5daba274dc83f5ede4e2acb67e2572693e564b3899bb2d69a60ee6d921c1af51f6834b792951808dfcba5f8841
    HEAD_REF develop
)

set(_is_shared TRUE)
if(VCPKG_LIBRARY_LINKAGE STREQUAL "static")
    set(_is_shared FALSE)
endif()

if("openmp" IN_LIST FEATURES)
    set(_use_)
endif()

vcpkg_check_features(OUT_FEATURE_OPTIONS FEATURE_OPTIONS
    FEATURES
        openmp       ENABLE_OPENMP
)



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
        ${FEATURE_OPTIONS}
)

vcpkg_install_cmake()
vcpkg_fixup_cmake_targets(CONFIG_PATH lib/cmake/camp
                          TARGET_PATH share/camp)
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
    # Note: Not tested
    file(REMOVE_RECURSE ${CURRENT_PACKAGES_DIR}/bin ${CURRENT_PACKAGES_DIR}/debug/bin)
else()
    set(_config_dir "${CURRENT_PACKAGES_DIR}/share/camp")

    # Move dll files from lib to bin directory
    file(MAKE_DIRECTORY ${CURRENT_PACKAGES_DIR}/bin )
    file(MAKE_DIRECTORY ${CURRENT_PACKAGES_DIR}/debug/bin )

    file(RENAME ${CURRENT_PACKAGES_DIR}/lib/camp.dll
                ${CURRENT_PACKAGES_DIR}/bin/camp.dll)

    file(RENAME ${CURRENT_PACKAGES_DIR}/debug/lib/camp.dll
                ${CURRENT_PACKAGES_DIR}/debug/bin/camp.dll)

    # Update paths to dlls in CMake config files
    foreach(_c  debug release)
        set(_f ${_config_dir}/campTargets-${_c}.cmake)
        file(READ ${_f} _fdata)
        string(REPLACE "lib/camp.dll" "bin/camp.dll" _fdata "${_fdata}")
        file(WRITE  ${_f} "${_fdata}")
    endforeach()
endif()

# Put the license file where vcpkg expects it
# Note: LICENSE occasionally cannot be found
if(EXISTS "${SOURCE_PATH}/LICENSE")
    file(INSTALL     ${SOURCE_PATH}/LICENSE 
         DESTINATION ${CURRENT_PACKAGES_DIR}/share/camp 
         RENAME      copyright)
endif()
