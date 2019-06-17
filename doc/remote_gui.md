# BioLinux Remote X access

BioLinux supports both VNC and freenx GUI X-windows access to a remote
VM. And you can use X programs through ssh, naturally.

## VNC

VNC is a ubiquitous remote access tool - always there, and easy to install/use.

In a nutshell:

Make sure vnc4server is installed on the VM.  Enable ports 5900, 5901 and 5902
on the VM. Run the server

        vnc4server -depth 24
           (set password)

Run the client on your desktop

       vncviewer -FullColor=1 HostIP:1

Where HostIP is the reachable host IP address or DNS name. Next, it 
may be necessary to start an X desktop, such as LXDE:

       startlxde

### Amazon EC2 ports

Create a security group for your instance that allows at least ports 
22,5900,5901 and 5902.

### Vagrant ports

You may need to add port forwarding to vagrant - as the testing system
does. I.e. add to the Vagrantfile:

       config.vm.forward_port('vnc0', 5900, 5900)
       config.vm.forward_port('vnc1', 5901, 5901)
       config.vm.forward_port('vnc2', 5902, 5902)

This is for testing, mostly. You do not need VNC on Vagrant/VirtualBox. Fire up
the GUI directly!

### VNC Security

Please note that VNC is not very secure - it has no proper key protection. You
can tunnel over ssh for improved security. Or use freenx instead.

## FreeNX

FreeNX is a fast version of the X protocol.

Make sure freenx is installed on the VM. CloudBioLinux comes with scripts
for setting up freenx.

(to be filled in)

## X over ssh

Normally you have ssh access to the remote VM. You can use X-windows programs
remotely, provided you have a local X server (always on Linux and OSX). Just
login with the -X switch

      ssh -X user@$hostIP

in the terminal type an X program, e.g.

      firefox

and the program should display locally (running remotely).

