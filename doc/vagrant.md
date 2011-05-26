# BioLinux and Vagrant

This document gives some additional information on using Vagrant with BioLinux.
[Vagrant][v1] is a convenient command line manager for VirtualBox. In conjunction
with the BioLinux fabric environment, any VirtualBox VM can be bootstrapped.

# Rolling your own

## Start from a BioLinux box

See the main README file for firing up a pre-installed BioLinux box.

## Start from scratch

Start with a standard downloadable prepared Vagrant box. For example a Debian
32-bits box prepared for Vagrant, or create one from scratch as is explained on
the [Vagrant web site][v1].

Next add the box to vagrant using a URL, or a file:

  vagrant box add debian_squeeze_32 debian_squeeze_32.box

and create your own version

  mkdir mynode
  cd mynode

Creates a ./Vagrantfile describing the VM.

  vagrant init debian_squeeze_32

Start the VM (which gets copied the first time):

  vagrant up

and login

  vagrant ssh

make sure you have enough disk space for the dir ~/VirtualBox\ VMs and
~/.vagrant, as this is where VMs are copied.

At this point a bare VM is running that will accept BioLinux installations. The next 
step is to pull the BioLinux tree and to run fab using the 'vagrant' host:

  fab -H vagrant -f /path/to/cloudbiolinux/fabfile.py install_biolinux

which uses the information from the local ./Vagrantfile.

[v1]: http://vagrantup.com/docs/base_boxes.html

  


