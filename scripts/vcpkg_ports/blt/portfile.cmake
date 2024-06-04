vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO llnl/blt
    REF v0.6.2
    SHA512 ca38639dc2bcbd9814e11257ce8d4c8b61e082dcc86ba6ceb339a6e473191b2008158014f06325a1a15b71897d5d42b0f16261dcdd6164081e415b11a06f3ff4
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
