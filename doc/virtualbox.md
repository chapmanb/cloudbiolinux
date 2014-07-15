# CloudBioLinux, VirtualBox and Vagrant

This document gives some additional information on using Vagrant with BioLinux.
[Vagrant][v1] is a convenient command line manager for VirtualBox. In conjunction
with the BioLinux fabric environment, any VirtualBox VM can be bootstrapped.

Note the current version of vagrant needs at least VirtualBox version 4.1.x.

## VirtualBox with vagrant

Add a base image to vagrant, and boot it up; community Vagrant boxes are available from
[http://vagrantbox.es][v3] and [BioLinux flavors][v4]:

        vagrant box add box_name http://path_to_the_image.box
        mkdir tmp/biolinux
        cd tmp/biolinux
        vagrant init box_name
        vagrant up

Run the fabfile, building CloudBioLinux:

        fab -H vagrant -f /path/to/cloudbiolinux/fabfile.py install_biolinux

Then build the box, renaming package.box to `cloudbiolinux_date` and
move it to a public webserver, such as Amazon S3:

        vagrant package
        mv package.box biolinux_20110122.box
        s3cmd put --acl-public --guess-mime-type biolinux_20110122.box
              s3://chapmanb/biolinux_20110122.box

[v3]: http://vagrantbox.es/
[v4]: http://biobeat.org/bionode

# Rolling your own

## Start from a BioLinux box

See the main README file for firing up a pre-installed BioLinux box.

## Start from scratch

Despite the extra work, starting from scratch may have advantages. For
one you have more control of the base install. Say for a different
version of Linux, a BSD kernel, or for install less software (do you
really need X/KDE/Gnome?), so you do not end up with an 8 GB VM, or for more
software and/or data pre-intstalled on a VM. 

The BioLinux setup is designed to be modular, to support multiple 
flavors (see the main README for an explanation of terms).

Start with a standard downloadable prepared Vagrant box. For example a Debian
32-bits box prepared for Vagrant, or create one from scratch as is explained on
the [Vagrant web site][v1].

Next add the virtualbox to vagrant using a URL, or box file:

          vagrant box add debian_squeeze_32 debian_squeeze_32.box

(boxes are available form [http://vagrantbox.es][v3] and
[http://biobeat.org/bionode][BioLinux flavors]) and create your own version

          mkdir myflavor
          cd myflavor

Creates a ./Vagrantfile describing the VM.

          vagrant init debian_squeeze_32

Have a look inside the Vagrantfile. The default should be fine now.

Start the VM (which gets copied the first time, which may take a while):

          vagrant up

and login

          vagrant ssh  # no password needed

make sure you have enough disk space (twice the box size) for the dir
~/VirtualBox\ VMs and ~/.vagrant, as this is where VMs are copied from the
original box file.

At this point a bare VM is running that will accept BioLinux installations. The
next step is to pull the BioLinux tree on your local system, and to run fab using the
vagrant host, using a minimal install target. E.g.

          export source=/path/to/cloudbiolinux

and

          fab -f $source/fabfile.py -H vagrant  -c $source/contrib/minimal/fabricrc_debian.txt install_biolinux:packagelist=$source/contrib/minimal/main.yaml

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
be fast.  Minimal runs the second time in under 20 seconds on a
basic laptop, as we do with a 'Minimal' install:

         ./test/test_vagrant --continue

For completeness, after a minimal install you can still install a full BioLinux
execute

        fab -H vagrant -f $source/fabfile.py install_biolinux

Once you have a working Virtual Box VM with vagrant, you can package it with

        vagrant package

and make the resulting .box file available for others to use.

Read the README for further information.

[v1]: http://vagrantup.com/docs/base_boxes.html

## Trouble shooting

### Guest additions

You may see an error

  [default] The guest additions on this VM do not match the install version of
  VirtualBox! This may cause things such as forwarded ports, shared
  folders, and more to not work properly. If any of those things fail on
  this machine, please update the guest additions and repackage the
  box.

  Guest Additions Version: 4.0.4
  VirtualBox Version: 4.1.0

this error may actually be caused by the Vbox Linux kernel drivers not having
been loaded! Fix

       modprobe vboxdrv

# Converting Vagrant images to VirtualBox and Eucalyptus images

(protocol steps tested in Ubuntu Natty)

## software pre-requisite

    sudo gem install vagrant
    sudo apt-get install cloud-utils

## Importing cloud biolinux VM to your system

    vagrant box add base 
    https://s3.amazonaws.com/cloudbiolinux/cbl_ubuntu_11_4_32_20110628.box
    vagrant init base
    vagrant up

## adding some missing components to the vagrant VMs 

    vagrant ssh
    sudo apt-get install gdm cloud-utils openssh
    sudo useradd -d /home/ubuntu -m ubuntusudo passwd ubuntu
    sudo shutdown -r now

in the graphical login after reboot get in with user:ubuntu / pass:ubuntu
go to System--->Administration--->Login Window to enable autologin

## VirtualBox Appliance

Virtual Appliances are pre-assemblied VM images configured for various purposes.

Open the Virtualbox GUI, you should see the VM added by vagrant - you can 
rename it to "Cloud BioLinux 32"

  File->Export Appliance

and distribute the .ova.

Anyone in any OS running Virtualbox can import the .ova with File->Import
Appliance.

# Making a Eucalyptus image from VirtualBox

Start with the Cloud BioLinux Virtualbox .vmdk (its location is in the VM
properties from the Virtualbox GUI). Resize the vmdk, since the size may be
40G, and the Eucalyptus image will have that size. 

Accordint to http://mtnbike.org/blog/?p=29 and the same here:
http://www.my-guides.net/en/content/view/122/26/

convert to raw .img 

    qemu-img convert -O raw CloudBioLinux-32bit-disk1.vmdk 
    CloudBioLinux-32bit-disk1.img

deploy to Eucalyptus via

    uec-publish-img CloudBioLinux-32bit-disk1.img

# VirtualBox, KVM or XEN?

There are more ways than one to virtualize machines on Linux.

Despite the attractions of vagrant and Virtualbox, as displayed here, we note
that Linux KVM may be a better choice for virtualization and testing of
CloudBioLinux, as Linux distributions support KVM out of the box, and KVM has
more Unix-like control.  See also the information for using KVM in
./doc/linux_kvm.md. 

For production environments check out XEN virtualization (XEN runs Amazon EC2).


