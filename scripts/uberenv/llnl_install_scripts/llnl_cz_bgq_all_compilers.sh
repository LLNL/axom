#!/bin/bash
#
# chang28 02-03-2017, script to set up BG/Q TPL build, for now it only works on Cyrus' branch
# chang28 02-07-2017, Cyrus is merging this branch to develop, no need to checkout his branch in the script anymore
# usage: 
#  cd {to directory with this script}
#   ./llnl_cz_bgq_all_compilers.sh

date
. /usr/local/tools/dotkit/init.sh
use python-2.7.3
use cmake-3.1
python ./llnl_bgq_uberenv_install_bgqos_all_compilers.py
date

