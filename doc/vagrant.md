# BioLinux and Vagrant

This document gives some additional information on using Vagrant with BioLinux.
[Vagrant][v1] is a convenient command line manager for VirtualBox. In conjunction
with the BioLinux fabric environment, any VirtualBox VM can be bootstrapped.

# Rolling your own

## Start from a BioLinux box

See the main README file for firing up a pre-installed BioLinux box.

## Start from scratch

Despite the extra work, starting from scratch may have some advantages. For one
you have more control of the base install. Say for a different version of
Linux, a BSD kernel, or for install less software (do you really need X?), so
you do not end up with an 8 GB VM. The BioLinux setup is designed to be
modular, to support multiple editions and flavors.

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

Start the VM (which gets copied the first time, which may take a while):

          vagrant up

and login

          vagrant ssh

make sure you have enough disk space for the dir ~/VirtualBox\ VMs and
~/.vagrant, as this is where VMs are copied.

At this point a bare VM is running that will accept BioLinux installations. The next 
step is to pull the BioLinux tree and to run fab using the vagrant
host, using a minimal install target. E.g.

          export source=/path/to/cloudbiolinux

and

          fab -f $source/fabfile.py -H vagrant  -c $source/contrib/minimal/fabricrc_debian.txt install_bare:packagelist=$source/contrib/minimal/main.yaml

which uses the information from the local ./Vagrantfile. 

The first time the minimal fabfile is run it updates the /etc/apt/sources (on
Debian-based systems), and a number of basic packages, including sudo, python,
chef. It may be the Linux kernel and support libraries get upgraded, if they
are in the dependency tree. Starting from a minimalistic Debian Vagrant box, the
BioLinux minimal install has an unpacked size under 1Gb. E.g.

        vagrant@vagrant-debian-squeeze:~$ df -h
        Filesystem            Size  Used Avail Use% Mounted on
        /dev/sda1              39G  804M   36G   3% /
        tmpfs                 188M     0  188M   0% /lib/init/rw
        udev                  184M  116K  184M   1% /dev
        tmpfs                 188M     0  188M   0% /dev/shm
        v-root                 51G   24G   27G  47% /vagrant

Despite the fact that running fabfile.py is destructive, i.e. it overwrites the
current install, it is reasonably safe as it mostly uses the underlying package
management system and dependency resolution. Rerunning a BioLinux fabfile can
be fast.  The minimal edition runs the second time in under 20 seconds on a
basic laptop.

For completeness, after a minimal install you can still install a full BioLinux
execute

        fab -H vagrant -f $source/fabfile.py install_biolinux

Read the README for further information.

[v1]: http://vagrantup.com/docs/base_boxes.html
