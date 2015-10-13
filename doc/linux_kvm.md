# CloudBioLinux and Linux KVM

This document gives some additional information on using BioLinux on Linux KVM.
KVM is a great virtualization environment, it is part of the Linux effort, will
work with the default Linux kernel. A running 64-bit Linux kernel can run both
32-bit and 64-bit VMs.  With most Linux distributions KVM comes out of the
box...

Together with the BioLinux fabric environment, any KMV VM can be
bootstrapped. Here we start from a bare Debian/Ubuntu
installation. With Debian you can install from a fresh download of
[stable](http://www.debian.org/releases/stable/) version burning it on
a CDROM.  The default install will do on, say, a 10GB partition - so
the rest can be used for LVM partitioning. Select ssh and the standard
system utilities.

# Install KVM

There are many web resources for installing KVM. On Debian derived systems:

      apt-get install kvm libvirt-bin virtinst bridge-utils \
                      qemu-kvm virt-manager libvirt-bin

and add your user to the kvm group. E.g.

    adduser username kvm

For non-Debian systems see, for example, [OpenSuse](http://doc.opensuse.org/documentation/html/openSUSE/opensuse-kvm/cha.kvm.requires.html).

Note that not all machines support virtualization, and if they do you
may need to switch it on in the BIOS (especially true on older
hardware). Check for CPU support with:

    egrep --color '(vmx|svm)' /proc/cpuinfo

Also note that some older Linux installations may need a kernel
upgrade. Start KVM with

    /etc/init.d/qemu-kvm start

or with later versions

    /etc/init.d/qemu-system-x86 start

# Create a bare VM

Download a live installation image file. For example fetch a standard
or network image from [[http://www.debian.org/][Debian]].

Reserve space on the disk partition - this should be enough for a
Debian install and updates.

      qemu-img create hda.img -opreallocation=metadata -ocluster_size=2M -f qcow2 10G

(settings suggested by Red Hat) and fire up the VM

      kvm debian-live-$(VER).img -hda hda.img -curses -no-reboot -serial pty

Alternatively use the netinstall. The CloudBioLinux integration test
system does something similar, starting from the smaller net install
of Debian Linux:

      qemu-system-x86_64 -enable-kvm -cdrom debian-$(VER)-amd64-netinst.iso -hda hda.img

hit ESC and optionally type 'install fb=false' to disable the frame buffer.
Fire up the installer. Note that the file system of the installer can be slow,
that speed is not representative for an installed VM later (with -enable-kvm).

With the base install, boot the new system

      qemu-system-x86_64 -enable-kvm -redir tcp:2222::22 -hda hda.img

and install ssh on the VM (it comes already on netinst)

      apt-get install openssh-server

so that ssh login works

      ssh -p 2222 biolinux@localhost

on user biolinux without a password (preferably using a key with
empty password). And give that user 'sudo bash'. This ssh and sudo
configuration is described in ./doc/private_cloud.md. After generating
the key 

      ssh -i ~/.ssh/biolinux -p 2222 biolinux@localhost

test run sudo without a password

      sudo bash

Not much to installing Linux with KVM! From this point onwards you can install
CloudBioLinux using the fabric file. Make a copy of the hda.img file,
so you can have the same starting point every time

      cp hda.img kvm_with_biolinux_login.img

The CloudBioLinux test script also starts from here.
Try the ./test/test_biolinux script to test drive the VM. test_biolinux
will install a CloudBioLinux flavor, and check whether the installation is
complete. Essentially with a running instance:

      ./test/test_biolinux -p 2222 -u biolinux -i ~/.ssh/biolinux 127.0.0.1

(note the use of 127.0.0.1 over localhost - this is because of a bug
in fabric - is this still true?).

# KVM tips

KVM is nice and powerful. It is used in many Cloud service providers.

For fast performance, it pays to install on a raw (LVM) partition,
get bridging sorted, and make sure hardware acceleration is in place.
Interesting goodies are the monitor (Crtl-Alt-2), virtsh, etc. See also
[http://www.linux-kvm.org/page/FAQ][kvm tips].

[kvm tips]: http://www.linux-kvm.org/page/FAQ
