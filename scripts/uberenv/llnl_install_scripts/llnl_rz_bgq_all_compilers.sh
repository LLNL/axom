#!/bin/bash
#
# chang28 03-17-2017, script to set up BG/Q TPL build on RZ
# usage: 
#  cd {to directory with this script}
#   ./llnl_rz_bgq_all_compilers.sh

date
. /usr/local/tools/dotkit/init.sh
use python-2.7.3
use cmake-3.1
python ./llnl_rz_bgq_uberenv_install_bgqos_all_compilers.py
date

