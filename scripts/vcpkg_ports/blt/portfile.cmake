vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO llnl/blt
    REF v0.5.1
    SHA512 b8971937d1f526d56afbb1ea664ca426840232ad99c85f0434b398f975719eb678c44e0f3dd24c9f36987403f685d4507adb2ae8a559c7bf28b79ac8de46a345
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
