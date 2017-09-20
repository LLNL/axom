#!/bin/bash
# axom_code_checkout.sh branch_name, if branch_name is NULL, develop is the default branch

echo  axom_code_checkout.sh 1.0
set -ev
 
rm -rf axom
BRANCH=${1:-develop}
echo $BRANCH
git clone ssh://git@cz-bitbucket.llnl.gov:7999/atk/axom.git --branch $BRANCH --recursive
cd axom
git branch
