# CloudBioLinux and Linux KVM

This document gives some additional information on using BioLinux on Linux KVM.

In conjunction with the BioLinux fabric environment, any KMV VM can be
bootstrapped.

# Install KVM

There are many web resources for installing KVM. On Debian derived systems:

      apt-get install kvm libvirt-bin virtinst bridge-utils \
                      qemu-kvm virt-manager libvirt-bin

and add your user to the kvm group. E.g.

      adduser user kvm
  
# Create a bare VM

Download a live installation image file. For example from

  http://cdimage.debian.org/cdimage/release/current-live/i386/usb-hdd/debian-live-$(VER)-i386-standard.img
  
Next reserve space on your disk partition (note, do not use qemu on an ext4 partition)

      qemu-img create hda.img -opreallocation=metadata -ocluster_size=2M -f qcow2 4G

(settings suggested by [http://docs.redhat.com/docs/en-US/Red_Hat_Enterprise_Linux/6/html/Technical_Notes/virt.html][redhat]) and fire up the VM

      kvm debian-live-$(VER)-i386-standard.img -hda hda.img -curses -no-reboot -serial pty

The CloudBioLinux integration test system does something similar, starting from
the smaller net install of Debian Linux:

      wget http://cdimage.debian.org/debian-cd/6.0.3/amd64/iso-cd/debian-6.0.3-amd64-netinst.iso
      qemu-system-x86_64 -cdrom debian-6.0.3-amd64-netinst.iso -hda hda.img

hit ESC and type 'install fb=false'. This will fire up the installer. With the
base install, boot the new system and set up ssh so it can be used without a
password (preferably using a key with empty password). From that point on you
can install CloudBioLinux using the fabric file.

This is also the image the test system fingerprints for further test installs. You can 
try the ./test/test_biolinux script to test drive the VM. test_biolinux will install a 
CloudBioLinux flavor, and check whether the installation is complete.


