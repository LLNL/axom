if(VCPKG_CMAKE_SYSTEM_NAME STREQUAL "WindowsStore")
    message(FATAL_ERROR "${PORT} does not currently support UWP")
endif()

include(vcpkg_common_functions)
vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO llnl/conduit
    REF v0.5.0
    SHA512 aae329cf7d0329b466e996f81695f4bee66e7732d0d7c49ffd00276ddee82326a261af3135ad1fc7b9903150cadbbaed17c7a1f25b4cc5352fdfed60ed7a7da1
    HEAD_REF master
    PATCHES
        "fix-setup-hdf5-vcpkg.patch"
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

# Move dll files to bin directory
foreach(_dll conduit conduit_blueprint conduit_relay)
    file(RENAME ${CURRENT_PACKAGES_DIR}/lib/${_dll}.dll
                ${CURRENT_PACKAGES_DIR}/bin/${_dll}.dll)
    file(RENAME ${CURRENT_PACKAGES_DIR}/debug/lib/${_dll}.dll
                ${CURRENT_PACKAGES_DIR}/debug/bin/${_dll}.dll)
endforeach()

# Update dll paths in config files from 'lib' to 'bin' directory
foreach(_build debug release)
    file(READ ${CURRENT_PACKAGES_DIR}/share/conduit/conduit-${_build}.cmake _conf_file)
    foreach(_dll conduit conduit_blueprint conduit_relay)
        string(REPLACE "\${_IMPORT_PREFIX}/debug/lib/${_dll}.dll"
                       "\${_IMPORT_PREFIX}/debug/bin/${_dll}.dll" _conf_file "${_conf_file}")
        string(REPLACE "\${_IMPORT_PREFIX}/lib/${_dll}.dll"
                       "\${_IMPORT_PREFIX}/bin/${_dll}.dll" _conf_file "${_conf_file}")
    endforeach()
    file(WRITE ${CURRENT_PACKAGES_DIR}/share/conduit/conduit-${_build}.cmake "${_conf_file}")
endforeach()

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


