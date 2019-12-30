if(VCPKG_CMAKE_SYSTEM_NAME STREQUAL "WindowsStore")
    message(FATAL_ERROR "${PORT} does not currently support UWP")
endif()

include(vcpkg_common_functions)
vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO llnl/blt
    REF develop
    SHA512 b4e21cd929773ac9ba596703a66e0b820dbe0778ef0626b5d39e9e3ac8a8b72347a4d5218da68df866c0635af1a2f027c92758937cf8890b7439374dace584f3
    HEAD_REF develop
)

file(MAKE_DIRECTORY
    ${CURRENT_PACKAGES_DIR}/include/
    ${CURRENT_PACKAGES_DIR}/share/blt
)

# Add an empty file to include to make vcpkg happy
file(WRITE ${CURRENT_PACKAGES_DIR}/include/empty.txt "")

# Copy files over to blt's 'share' directory
set(_src_dir ${SOURCE_PATH})
set(_dest_dir ${CURRENT_PACKAGES_DIR}/share/blt)

foreach(_f cmake LICENSE RELEASE-NOTES.md SetupBLT.cmake thirdparty_builtin tests )
  file(COPY ${_src_dir}/${_f} DESTINATION ${_dest_dir} )
endforeach()

# Rename the LICENSE file to 'copyright'
file(RENAME ${_dest_dir}/LICENSE ${_dest_dir}/copyright)
