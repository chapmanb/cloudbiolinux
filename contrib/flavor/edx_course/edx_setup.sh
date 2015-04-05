#!/bin/bash
set -eu -o pipefail

# Vagrant installation script for creating VM to use in the
# edX variant analysis course, PH525.6x:
# https://www.edx.org/course/case-study-variant-discovery-genotyping-harvardx-ph525-6x

# Base system
sudo apt-get update
sudo apt-get install -y build-essential zlib1g-dev wget curl python-setuptools git
sudo apt-get install -y openjdk-7-jdk openjdk-7-jre ruby libncurses5-dev libcurl4-openssl-dev libbz2-dev unzip pigz bsdmainutils
sudo apt-get install -y python-pip python-dev
sudo pip install fabric

# CloudBioLinux
cd /vagrant
sudo chown -R vagrant:vagrant /usr/local
git clone https://github.com/chapmanb/cloudbiolinux.git
fab -H localhost install_biolinux:flavor=edx_course
# Need to hack cmake (to build 32bit): brew uninstall --force cmake && brew install cmake --without-docs
#
# GEMINI
wget https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py
python gemini_install.py /usr/local /usr/local/share/gemini --nosudo --nodata --notools

# Copy snpEff hg19 file from existing bcbio installation

# cleanup
rm -rf ~/.cache/Homebrew
rm -rf /usr/local/share/gemini/gemini
/usr/local/share/gemini/anaconda/bin/conda clean --yes --tarballs
sudo apt-get clean
sudo rm -rf /var/lib/apt/lists/*

# On VM -- Compact space: http://andrewdeponte.com/2013/10/29/shrinking-vagrant-linux-boxes.html
sudo dd if=/dev/zero of=wipefile bs=1024x1024; rm wipefile
# Afterwards, clone and export:
# vagrant halt
# virtualbox
#  - Select vagrant image
#  - Machine -> Clone to new VirtualMachine
#  - File -> Export as OVA
