#!/bin/bash
########################################
# use uberenv to build third party libs.
########################################

USAGE="$0 [-c <compiler> -p <prefix>]\nCheck uberenv compiler.yaml file for list of valid compilers."

# Set default values.
OPT_COMPILER="gcc@4.9.3"
b_hasPrefix=0

while getopts ':c:p:' opt
do
  case $opt in
    c) OPT_COMPILER=$OPTARG;;
    p) OPT_PREFIX=$OPTARG
       b_hasPrefix=1;;
   \?) printf "Error: invalid option.\nUsage: $USAGE\n"
       exit 1;;
  esac
done

if [ $b_hasPrefix -ne 1 ]; then
  printf "Must specify a prefix with '-p'.\n"
  exit 1
fi

export LC_TPL_PATH=$OPT_PREFIX
# create new library directory
python ../uberenv.py --prefix $OPT_PREFIX --spec %$OPT_COMPILER

# Update permissions on installed files.
./lc_install_set_perms.sh
