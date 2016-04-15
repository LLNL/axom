#!/bin/bash
########################################
# use uberenv to build third party libs.
########################################

USAGE="$0 [-c <compiler> -p <prefix> -m <mirror>]\nCheck uberenv compiler.yaml file for list of valid compilers."

# Set default values.
OPT_COMPILER="gcc@4.9.3"
b_hasPrefix=0

# use LC mirror if it exists
OPT_MIRROR=/usr/gapps/asctoolkit/thirdparty_libs/mirror
if [ -d  $OPT_MIRROR ]; then
  b_hasMirror=1
else
  OPT_MIRROR=""
  b_hasMirror=0
fi

while getopts ':c:p:m:' opt
do
  case $opt in
    c) OPT_COMPILER=$OPTARG;;
    p) OPT_PREFIX=$OPTARG
       b_hasPrefix=1;;
    m) OPT_MIRROR=$OPTARG
       b_hasMirror=1;;
   \?) printf "Error: invalid option.\nUsage: $USAGE\n"
       exit 1;;
  esac
done

if [ $b_hasPrefix -ne 1 ]; then
  printf "Must specify a prefix with '-p'.\n"
  exit 1
fi

if [ $b_hasMirror ]; then
  FLAG_MIRROR="--mirror $OPT_MIRROR"
else
  FLAG_MIRROR=""
fi


export LC_TPL_PATH=$OPT_PREFIX
# create new library directory
python ../uberenv.py --prefix $OPT_PREFIX --spec %$OPT_COMPILER $FLAG_MIRROR

# Update permissions on installed files.
./lc_install_set_perms.sh
