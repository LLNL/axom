#
# Generate source
# slic/src/interface/CMakeLists.txt
#
add_shroud( YAML_INPUT_FILE slic_shroud.yaml
            C_FORTRAN_OUTPUT_DIR c_fortran
)

uncrustify_shroud()

