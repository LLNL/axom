#!/bin/bash
#
# chang28 02-03-2017, script to set up BG/Q TPL build, for now it only works on Cyrus' branch
# usage: 
#  cd {to directory with this script}
#   ./llnl_bgq_toss3_all_compilers.sh

date
. /usr/local/tools/dotkit/init.sh
use python-2.7.3
use cmake-3.1
git checkout task/cyrush/2016_12_update_conduit
python ./llnl_bgq_uberenv_install_bgqos_all_compilers.py
date

