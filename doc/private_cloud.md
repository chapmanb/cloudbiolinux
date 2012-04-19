# Private Cloud and CloudBioLinux

CloudBioLinux can be used to create a private Cloud for
Bioinformatics. Essentially, all you need is ssh access to a VM
running somewhere. This VM should be a clean install of Linux.  With
CloudBioLinux Debian and Ubuntu distributions are supported best.

## Start a VM

Start a VM and make sure there is a network defined, and ssh running.
On the VM

       dhclient -v
       ifconfig
       ps xa|grep ssh
       ssh localhost

Check the network (e.g. with Debian)

       apt-get update
       apt-get install vim

You should be able to use the IP address to login from your desktop

       ssh user@VM_IP_address

## Get password free ssh access

The CloudBioLinux fabric tools work best when you have password free
login. If you can login to the remote with

       ssh user@VM_IP_address

you are set. Otherwise, create a password free ssh key. To achieve
this, see the many Internet resources, e.g.
http://www.mtu.net/~engstrom/ssh-agent.php.

## Install sudo without password

Install the sudo program. Next, edit /etc/sudoers with the 'visudo'
command, and add the line

       user ALL=NOPASSWD: /bin/bash

where user is your VM user login name. Alternatively add user to the sudo
group.

Now try:

       sudo bash

and you should be root, without a password.

## Install CloudBioLinux

See the README for installing CloudBioLinux and fabric.

## Run fabric

Now you should be set! To install BioLinux

       fab -f $source/fabfile.py -H user@$VM_IP_address -c $fabricrc install_biolinux:packagelist=$packagelist

Where source points to the checked out source tree, e.g.

       export source=$HOME/izip/git/opensource/debian/biolinux

For example, to install the Minimal flavor on Debian stable on a VM
running on IP 192.168.64.105:

       fab -f $source/fabfile.py -H user@192.168.64.105 \
       -c $source/contrib/minimal/fabricrc_debian.txt \
       install_biolinux:packagelist=$source/contrib/minimal/main.yaml

CloudBioLinux shows the following output. First it sets up the
environment

        [192.168.64.105] Executing task 'install_biolinux'
        cloudbiolinux WARNING: Skipping fabricrc.txt as distribution is already defined
        cloudbiolinux DEBUG: Minimal Edition 1.0.1
        cloudbiolinux INFO: This is a minimal
        cloudbiolinux INFO: This is a Base Flavor - no overrides
        cloudbiolinux INFO: Distribution debian
        cloudbiolinux INFO: Debian setup
        cloudbiolinux DEBUG: Debian-shared setup
        cloudbiolinux DEBUG: Source=squeeze
        cloudbiolinux DEBUG: Checking target distribution debian
        [192.168.64.105] run: cat /proc/version
        [192.168.64.105] out: Linux version 2.6.32-5-amd64 (Debian 2.6.32-31) (ben@decadent.org.uk) (gcc version 4.3.5 (Debian 4.3.5-4) ) #1 SMP Mon Mar 7 21:35:22 UTC [192.168.64.105] out:
        [192.168.64.105] out:
        cloudbiolinux INFO: Now, testing connection to host...
        cloudbiolinux INFO: Connection to host appears to work!
        cloudbiolinux DEBUG: Expand paths
        cloudbiolinux INFO: packagelist=/home/user/izip/git/opensource/debian/biolinux/contrib/minimal/main.yaml
        cloudbiolinux INFO: Meta-package information
        cloudbiolinux INFO: minimal,ruby
        cloudbiolinux INFO:
        cloudbiolinux INFO: Target=None

Here it modifies the source file for apt-get, as well as keys:

        cloudbiolinux DEBUG: _setup_apt_sources /etc/apt/sources.list.d/cloudbiolinux.list Minimal Edition
        [192.168.64.105] sudo: touch /etc/apt/sources.list.d/cloudbiolinux.list
        [192.168.64.105] sudo: echo '# This file was modified for Minimal Edition' >> /etc/apt/sources.list.d/cloudbiolinux.list
        cloudbiolinux DEBUG: Source deb http://ftp.nl.debian.org/debian/ squeeze main contrib non-free
        [192.168.64.105] sudo: echo 'deb http://ftp.nl.debian.org/debian/ squeeze main contrib non-free' >> /etc/apt/sources.list.d/cloudbiolinux.list
        cloudbiolinux DEBUG: Source deb http://ftp.nl.debian.org/debian/ squeeze-updates main contrib non-free
        [192.168.64.105] sudo: echo 'deb http://ftp.nl.debian.org/debian/ squeeze-updates main contrib non-free' >> /etc/apt/sources.list.d/cloudbiolinux.list
        [192.168.64.105] sudo:
        cloudbiolinux INFO: Update GPG keys for repositories
        cloudbiolinux INFO: Update and install all packages
        [192.168.64.105] sudo: apt-get update
        [192.168.64.105] out: Hit http://ftp.nl.debian.org squeeze Release.gpg

and starts installing packages

        cloudbiolinux INFO: Updating 26 packages
        [192.168.64.105] sudo: apt-get -y --force-yes install ruby1.8 ruby1.8-dev ruby1.9.1 ruby1.9.1-dev axel less openssh-server rsync screen sudo tar unzip bzr cvs darcs git-core mercurial subversion vim cmake g++ gcc gfortran make patch swig
        [192.168.64.105] out: Reading package lists... Done
        [192.168.64.105] out: Building dependency tree
        [192.168.64.105] out: Reading state information... Done
        [192.168.64.105] out: gcc is already the newest version.
        [192.168.64.105] out: gcc set to manually installed.
        [192.168.64.105] out: less is already the newest version.
        [192.168.64.105] out: less set to manually installed.
        [192.168.64.105] out: make is already the newest version.
        (etc, etc)

Finally some clean ups

        [192.168.64.105] sudo: apt-get clean
        cloudbiolinux INFO: Target=unknown; Edition=Minimal Edition; Flavor=Base Flavor - no overrides

write an entry in the log file

        [192.168.64.105] sudo: date +"%D %T - Updated Target=unknown; Edition=Minimal Edition; Flavor=Base Flavor - no overrides" >> /var/log/biolinux.log
        [192.168.64.105] run: uname -m
        [192.168.64.105] out: x86_64
        [192.168.64.105] out:
        cloudbiolinux INFO: Reading /home/user/izip/git/opensource/debian/biolinux/config/custom.yaml
        cloudbiolinux DEBUG: Packages:
        cloudbiolinux DEBUG:
        cloudbiolinux INFO: Cleaning up space from package builds
        [192.168.64.105] sudo: rm -rf .cpanm
        [192.168.64.105] sudo: rm -f /var/crash/*

And it is done. Minimal has no post-installation configuration, but
that is easy to add.

