#! /bin/sh 
#
# Install the biolinux-phylogeny on a host
#
# Usage:
#
#   ./contrib/flavor/phylogeny/install_host.sh user@hostname
#
# where
#
#   user:      the password-less login name (see doc/private_cloud.md)
#   hostname:  the hostname, or IP
#
#
# Copyright (C) 2011 Pjotr Prins <pjotr.prins@thebird.nl>

host=$1
if [ ! -e fabfile.py -o -z $host ]; then
  echo Usage:
  echo
  echo   ./contrib/flavor/phylogeny/install_host.sh user@hostname
fi
source=`pwd`
fabricrc=$source/contrib/flavor/phylogeny/fabricrc_debian.txt
packagelist=$source/contrib/flavor/phylogeny/main.yaml
fab -f $source/fabfile.py -H $host -c $fabricrc install_biolinux:packagelist=$packagelist
