#!/usr/bin/env bash
#------------------------------------------------------------------------------
# This script installs client-side hooks
#------------------------------------------------------------------------------
basedir=`git rev-parse --show-toplevel`
hooksdir="$basedir/.git/hooks/"
cp -v $basedir/scripts/gitsetup/hooks/commit-msg $hooksdir
