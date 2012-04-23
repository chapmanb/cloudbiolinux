# CloudBioLinux and Linux KVM

This document gives some additional information on using BioLinux on Linux KVM.
KVM is a nice virtualization environment, it is part of the Linux effort, will
work with the default Linux kernel, and comes with brilliant tools. A 64-bit
Linux can run both 32-bit and 64-bit VMs.

Linux VMs are useful - especially when you deal with software upgrade paths. A VM 
can act like a try-before-you-buy environment. But it is more that that - you
can move VMs from one machine to another, share machines with multiple
installed systems, etc. If you can install Linux for a server, it is probably
a good idea to use KVM. And it just comes out of the box...

In conjunction with the BioLinux fabric environment, any KMV VM can be
bootstrapped. Here we start from a basic Debian/Ubuntu/Mint installation. With
Debian you can install from a fresh download of
[stable](http://www.debian.org/releases/stable/) version burning it on a CDROM.
The default install will do on, say, a 6GB partition - so the rest can be used
for LVM partitioning. Select ssh and the standard system utilities.

# Install KVM

There are many web resources for installing KVM. On Debian derived systems:

      apt-get install kvm libvirt-bin virtinst bridge-utils \
                      qemu-kvm virt-manager libvirt-bin

and add your user to the kvm group. E.g.

      adduser user kvm

For non-Debian systems see, for example, [OpenSuse](http://doc.opensuse.org/documentation/html/openSUSE/opensuse-kvm/cha.kvm.requires.html).

Note that not all machines support virtualization, and if they do you may need
to switch it on in the BIOS (especially true on older hardware). Check for CPU support with:

  egrep --color '(vmx|svm)' /proc/cpuinfo

Also note that on older Linux installations may need a kernel upgrade.

# Create a bare VM

Download a live installation image file. For example from

  http://cdimage.debian.org/cdimage/release/current-live/i386/usb-hdd/debian-live-$(VER)-i386-standard.img

Next reserve space on your disk partition (note, do not use qemu on an ext4 partition)

      qemu-img create hda.img -opreallocation=metadata -ocluster_size=2M -f qcow2 4G

(settings suggested by [redhat][redhat]) and fire up the VM

      kvm debian-live-$(VER)-i386-standard.img -hda hda.img -curses -no-reboot -serial pty

The CloudBioLinux integration test system does something similar, starting from
the smaller net install of Debian Linux:

      wget http://cdimage.debian.org/debian-cd/6.0.3/amd64/iso-cd/debian-6.0.3-amd64-netinst.iso
      qemu-system-x86_64 -cdrom debian-6.0.3-amd64-netinst.iso -hda hda.img

hit ESC and type 'install fb=false'. This will fire up the installer. With the
base install, boot the new system

      qemu-system-x86_64 -enable-kvm -redir tcp:2222::22 -hda hda.img

and set up ssh on the VM

      apt-get install openssh-server

so this works

      ssh -p 2222 biolinux@localhost

so it can be used on a user without a password (preferably using a key with
empty password). Also give that use 'sudo bash'. This ssh and sudo
configuration is described in ./doc/private_cloud.md.

From that point onwards you can install CloudBioLinux using the fabric file.

This is also the image the test system fingerprints for further test installs.

You can try the ./test/test_biolinux script to test drive the VM. test_biolinux
will install a CloudBioLinux flavor, and check whether the installation is
complete. Essentially with a running instance

      ./test/test_biolinux -p 2222 -u biolinux 127.0.0.1

(note the use of 127.0.0.1 over localhost - this is because of a bug
in fabric).

# KVM tips

KVM is powerful. Performance-wise it pays to install on a raw (LVM) partition,
get bridging sorted, and make sure hardware acceleration is in place.
Interesting goodies are the monitor (Crtl-Alt-2), virtsh, etc. See also
[http://www.linux-kvm.org/page/FAQ][kvm tips].

[kvm tips]: http://www.linux-kvm.org/page/FAQ
[redhat]: http://docs.redhat.com/docs/en-US/Red_Hat_Enterprise_Linux/6/html/Technical_Notes/virt.html
