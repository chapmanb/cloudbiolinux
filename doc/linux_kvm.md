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

Download a CDROM installation file. For example from

  http://cdimage.debian.org/cdimage/release/current-live/i386/usb-hdd/debian-live-$(VER)-i386-standard.img
  
Next reserve space on your disk partition

      qemu-img create hda.img -f qcow 4G

and fire up the installer

      kvm debian-live-$(VER)-i386-standard.img -hda hda.img -curses -no-reboot -serial pty
  



