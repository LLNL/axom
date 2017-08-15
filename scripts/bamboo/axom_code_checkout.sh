#!/bin/bash
# axom_code_checkout.sh branch_name

echo  axom_code_checkout.sh 1.0
set -ev
 
rm -rf axom
 
git clone ssh://git@cz-bitbucket.llnl.gov:7999/atk/axom.git --branch $1 --recursive
cd axom
git branch
