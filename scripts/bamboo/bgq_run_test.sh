echo "compiler_flag_unused_var"
srun -n 4 -ppsmall ../../atk_build/tests/compiler_flag_unused_var 
wait $!
echo "compiler_flag_omp_pragma"
srun -n 4 -ppsmall ../../atk_build/tests/compiler_flag_omp_pragma 
wait $!
echo "compiler_flag_unused_param"
srun -n 4 -ppsmall ../../atk_build/tests/compiler_flag_unused_param
wait $!
echo "compiler_flag_strict_aliasing"
srun -n 4 -ppsmall ../../atk_build/tests/compiler_flag_strict_aliasing
wait $!
echo "blt_mpi_smoke"
srun -n 4 -ppsmall ../../atk_build/tests/blt_mpi_smoke
wait $!
echo "compiler_flag_uninitialized"
srun -n 4 -ppsmall ../../atk_build/tests/compiler_flag_uninitialized
wait $!
echo "lumberjack_speedTest"
srun -n 4 -ppsmall ../../atk_build/tests/lumberjack_speedTest
wait $!
echo "blt_gtest_smoke"
srun -n 4 -ppsmall ../../atk_build/tests/blt_gtest_smoke
wait $!
echo "common_fileUtilities"
srun -n 4 -ppsmall ../../atk_build/tests/common_fileUtilities
wait $!
echo "common_types"
srun -n 4 -ppsmall ../../atk_build/tests/common_types
wait $!
echo "lumberjack_RootCommunicator"
srun -n 4 -ppsmall ../../atk_build/tests/lumberjack_RootCommunicator
wait $!
echo "lumberjack_TextEqualityCombiner"
srun -n 4 -ppsmall ../../atk_build/tests/lumberjack_TextEqualityCombiner
wait $!
echo "common_TaskTimer"
srun -n 4 -ppsmall ../../atk_build/tests/common_TaskTimer
wait $!
echo "conduit_smoke"
srun -n 4 -ppsmall ../../atk_build/tests/conduit_smoke
wait $!
echo "slic_asserts"
srun -n 4 -ppsmall ../../atk_build/tests/slic_asserts
wait $!
echo "slam_AccessingRelationDataInMap"
srun -n 4 -ppsmall ../../atk_build/tests/slam_AccessingRelationDataInMap
wait $!
echo "lumberjack_Lumberjack"
srun -n 4 -ppsmall ../../atk_build/tests/lumberjack_Lumberjack
wait $!
echo "common_config"
srun -n 4 -ppsmall ../../atk_build/tests/common_config
wait $!
echo "lumberjack_Message"
srun -n 4 -ppsmall ../../atk_build/tests/lumberjack_Message
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/primal_bvhtree
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/fmt_smoke
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/slic_fmt
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/primal_boundingbox
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/primal_intersection
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/primal_orient
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/primal_intersection_impl
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/primal_numeric_array
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/primal_squared_distance
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/primal_point
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/primal_vector
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/primal_triangle
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/slam_DynamicVariableRelation
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/slam_IndirectionSet
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/slam_NullSet
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/slam_ModularInt
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/slam_Map
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/slam_utilities
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/slam_RangeSet
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/slam_Set
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/sidre_buffer
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/slam_StaticConstantRelation
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/sidre_buffer_unit
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/slam_WindowedRangeSet
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/slam_StaticVariableRelation
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/sidre_class
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/sidre_datastore
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/sidre_external
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/sidre_opaque
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/sidre_native_layout
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/sidre_datastore_unit
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/mint_mesh_extent
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/sidre_smoke
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/sidre_group
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/sidre_view
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/spio_externalWriteRead
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/spio_irregularWriteRead
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/spio_parallelWriteRead
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/spio_basicWriteRead
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/slam_tinyHydro_unitTests
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/quest_octree
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/quest_signed_distance
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/quest_spatial_octree
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/quest_morton
wait $!
srun -n 4 -ppsmall ../../atk_build/tests/quest_inout_octree
