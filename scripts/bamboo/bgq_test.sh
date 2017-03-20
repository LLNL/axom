# chang28 01-23-2017, this is the script that kicks off test
# chang28 01-25-2017, 
# chang28 01-26-2017, 
# chang28, 01-31-2017, change the script to generate bgq_autogen_test.sh on the fly to run this script please cd asctoolkit and ./scripts/bamboo/bgq_test.sh
# chang28, 02-28-2017, project_test_commands containsa a list of all test names and commands for the configuration. The grep command reads each line of project_test_commands.txt and produces two lines. 

echo "running bgq_test.sh version 1.00"
echo "cd atk_build"
cd atk_build
#grep -v ^# project_test_commands.txt | awk '{$1 = ""; printf "echo "; print; printf "srun -N1 -n1 -ppsmall " ;print}' > bgq_autogen_test.sh
find ./tests -print > bgq_test_lists.txt
grep -v ^# bgq_test_lists.txt | awk '{printf "echo "; print; printf "srun -N1 -n1 -ppsmall " ;print}' > bgq_autogen_test.sh
chmod a+x bgq_autogen_test.sh
