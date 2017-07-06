#
# Generate source
# spio/src/interface/CMakeLists.txt
#
set(YAML_INPUT_FILE spio_shroud.yaml)

add_shroud(
    YAML_INPUT_FILE      ${YAML_INPUT_FILE}
    YAML_OUTPUT_DIR      yaml
    C_FORTRAN_OUTPUT_DIR c_fortran
)
