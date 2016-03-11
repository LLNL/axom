#!/bin/bash
export LC_TPL_PATH=/usr/gapps/asctoolkit/thirdparty_libs/`date +%Y_%m_%d`
if [ ! -d "$LC_TPL_PATH" ]; then
  echo "Creating tpl path $LC_TPL_PATH"
  mkdir $LC_TPL_PATH
fi
