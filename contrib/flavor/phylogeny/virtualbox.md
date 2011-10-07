= Running the Phylogeny VM in Virtualbox

== Install Virtualbox

Download and install Virtualbox from https://www.virtualbox.org/

Start Virtualbox.

== Download and install the BioLinux image

Fetch a BioLinux for virtualbox using the provided link, e.g.

        wget http://hostname/biolinux-phylogeny-virtualbox.vmdk

Add this file to VirtualBox by selecting 'New', choose a name and
select Debian Linux. After setting RAM, select 'Use existing hard disk' and
select the downloaded .vmdk file.

== Start up the image

Select the image and press 'Start'.

== Using the VM

After starting the VM you get a desktop. The user name is 'vagrant', as well
as the password. You can get root with the same password. Open a terminal
with LXterminal (from the menu). Run tools from the command line, e.g.

       raxmlHPC
       mb-mpi

you can have root, so installing new software is possible using apt-get.

== Modifying the virtual hardware

You can configure your image to run on multiple CPUs. Right-click on the image
icon in virtualbox and 'Settings'. Configure the hardware.

== Use ssh/scp

Openssh-server should be running on the VM

        apt-get install openssh-server

The VM is visible on a local network inside your PC. The VM network card is
eth1. Start a terminal as root and initialize

        dhcpclient -v eth1
        /sbin/ifconfig

from your PC use ssh to the listed IP address. 

When DHCP fails, reconfigure networking in VirtualBox. Right-click on the image
icon of the VM, click 'Settings', and configure the network to use 'bridge'.
And again run DHCP.

== Copying files to the image

With scp files can be copied to the image. An alternative approach, since the
VM has network, is to provide files through a webserver, or Amazon S3. E.g.

         wget http://webserver/file.tgz

