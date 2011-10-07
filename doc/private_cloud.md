= Private Cloud and CloudBioLinux

CloudBioLinux can be used to create a private Cloud for
Bioinformatics. Essentially, all you need is ssh access to a VM
running somewhere. This VM should be a clean install of Linux.  With
CloudBioLinux Debian and Ubuntu distributions are supported best.

== Start a VM

Start a VM and make sure there is a network defined, and ssh running.
On the VM

       /sbin/ifconfig
       ps xa|grep ssh
       ssh localhost

You should be able to use the IP address to login from your desktop

       ssh user@VM_IP_address

== Get password free ssh access

The CloudBioLinux fabric tools work best when you have password free
login. If you can login to the remote with 

       ssh user@VM_IP_address

you are set. Otherwise, create a password free ssh key. To achieve
this, see the many Internet resources, e.g.
http://www.mtu.net/~engstrom/ssh-agent.php.

== Install sudo without password

Install the sudo program. Next, edit /etc/sudoers with the 'visudo' 
command, and add the line

       user ALL=NOPASSWD: /bin/bash

where user is your Fabric username. Alternatively add user to the sudo
group.

Now try:

       sudo bash

and you should be root, without a password.
