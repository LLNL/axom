# chang28 01-23-2017, this is the script that kicks of test
echo "running bgq_test.sh"
echo "srun -N1 -n1 -ppsmall ../../atk_build/tests/slic_fmt"
srun -N1 -n1 -ppsmall ../../atk_build/tests/slic_fmt
