vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO llnl/blt
    REF v0.5.2
    SHA512 63e06282483985ae056e4a1557f249a9629130a4f5826c33ef3dbb7b8b1fc1760898fa89abd9734c3ab740aaf253e7284bad6aa342b92286ece810afe62350c2
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
