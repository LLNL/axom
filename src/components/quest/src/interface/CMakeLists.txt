#
# Generate source
# quest/src/interface/CMakeLists.txt
#
set(YAML_INPUT_FILE quest_shroud.yaml)

add_shroud(
    YAML_INPUT_FILE      ${YAML_INPUT_FILE}
    C_FORTRAN_OUTPUT_DIR c_fortran
)
install(FILES ${CMAKE_Fortran_MODULE_DIRECTORY}/axom_quest.mod DESTINATION lib/fortran)
