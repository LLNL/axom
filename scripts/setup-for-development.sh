#!/usr/bin/env bash

##
## Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
##
## Produced at the Lawrence Livermore National Laboratory.
##
## LLNL-CODE-741217
##
## All rights reserved.
##
## This file is part of Axom.
##
## For details about use and distribution, please read axom/LICENSE.
##

##------------------------------------------------------------------------------
## Script to setup dev environment.
##------------------------------------------------------------------------------

## get current path
basedir=`pwd`

reset=$(tput sgr0)
bold=$(tput bold)

echo "$bold Setting up development environment...$reset"
# ensure we are inside the repository
cd "${BASH_SOURCE%/*}/.."

GITSETUP="scripts/gitsetup"

current_directory=`pwd`
echo "$bold Current directory: $current_directory $reset"
echo "$bold == User Setup == $reset"
$GITSETUP/setup-git-user.sh &&
echo "$bold == Setup Git Editor == $reset" &&
$GITSETUP/setup-git-editor.sh &&
echo "$bold == Setup Git Aliases == $reset" &&
$GITSETUP/setup-git-aliases.sh &&
echo "$bold == Setup Git Colors == $reset" &&
$GITSETUP/setup-git-colors.sh &&
echo "$bold == Installing client-side Git hooks == $reset" &&
$GITSETUP/setup-git-hooks.sh &&
echo "$bold == Useful Tips/Suggestions == $reset" &&
$GITSETUP/tips.sh
