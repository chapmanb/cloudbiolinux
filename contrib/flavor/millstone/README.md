The Millstone flavor of cloudbiolinux eases deploying the Church Lab's
[Millstone](http://churchlab.github.io/millstone/) software,
a platform for genome design and analysis.

This document is intended for developers. Most users will want to clone an
existing, pre-packaged Amazon Machine Image (AMI), Those instructions are
available at the main Millstone webpage here: <http://churchlab.github.io/millstone/>.

The following instructions provided steps for deploying a (possibly) modified
version of Millstone to AWS. The instructions are specific to Amazon, but the
steps can be modified for development in different VM environments.

This document describes how to deploy an instance and create an AMI.
See the main cloudbiolinux README.md for general usage instructions for
cloudbiolinux.

### Usage

The following steps provide instructions for configuring a master, worker, or
hybrid instance of Millstone. The MASTER instance runs the main Django web
server as well as the Postgresql database and hosts the central RabbitMQ
messaging system for distributing computation. A WORKER instance runs a minimal
web server to allow the user to connect it to the master instance, and then
run asynchronous tasks as they become available.

When executing a fabric install, be sure to set either `MASTER=1` or `WORKER=1`, or both in environment variables to specify the role (Master or Worker) of the VM. Setting both runs the web server and background workers on the same machine. See below for example commands.

*NOTE*: The split master/worker implementation is currenlty under development. Only the combined setup has been tested.

#### Steps

##### Local requirements

On your local machine, install Python-Fabric (and its dep PyYAML), which are
used by cloudbiolinux to manage deploying to a remote machine.

    pip install fabric
    pip install pyyaml

##### Commands to load VMs

The following commands are run from the top-level cloudbiolinux directory
on your local machine and convert a clean EC2 instance to running a combined,
master, or worker configuration of Millstone, respectively.

*NOTE*: Remember to replace both `SSH_KEY_FILE` and `ec2-xx-xx-xx-xx.compute-1.amazonaws.com` with your ec2 private key and the correct instance url.

###### Master/Worker combined VM:

    MASTER=1 WORKER=1 fab -i SSH_KEY_FILE -u ubuntu -H ec2-xx-xx-xx-xx.compute-1.amazonaws.com install_biolinux:flavor=millstone

###### Master VM only (IN PROGRESS)

    MASTER=1 fab -i SSH_KEY_FILE -u ubuntu -H ec2-xx-xx-xx-xx.compute-1.amazonaws.com install_biolinux:flavor=millstone

###### Worker VM (IN PROGRESS)

    WORKER=1 fab -i SSH_KEY_FILE -u ubuntu -H ec2-xx-xx-xx-xx.compute-1.amazonaws.com install_biolinux:flavor=millstone

*DO NOT RESTART* if you want to package an Amazon Machine Image (AMI).
Restarting will bootstrap the software and data, resulting in potentially
insecure/malfunctional configurations. Immediately follow the steps in the
section on Building an AMI below.

##### Connecting Workers to Master

*NOTE*: Not yet supported. Draft:

If you configured separate master and worker instances, you will need to tell
the worker instances about the master instance.

The master info is located on the master instance at:

    http://ec2-xx-xx-xx-xx.compute-1.amazonaws.com/ec2/info

This shows its hostname and password, two important information that need to
put into the worker instance.

On worker instance, add the above information here:

    http://ec2-xx-xx-xx-xx.compute-1.amazonaws.com/worker

You need to update hostname and password from the master instance, so that the
worker processes can connect properly. No need to change other fields. That
page then allows you to test the connection.

The master will automatically put async tasks on its RabbitMQ and registered
workers will consume them. Workers will write back to the master's Postgresql
database.

##### Building an Amazon Machine Image (AMI)

First, clean up remaining SSH keys and history files.

    fab -c contrib/flavor/millstone/ec2.txt -i SSH_KEY_FILE -u ubuntu -H ec2-xx-xx-xx-xx.compute-1.amazonaws.com install_biolinux:flavor=millstone,target=clean

Go to EC2 Dashboard inside AWS Console <<a href="https://console.aws.amazon.com/">console.aws.amazon.com</a>>. Find the instance on which we just installed our software, and select "Create Image" in "Actions" tab.

##### Creating VM from AMI

Inside EC2 Dashboard, select create instance and choose AMI built in previous section. Initial bootstrap will take about five minutes for master VM and one minute for worker VM, before the VM becomes responsive.

If you instantiate separate master and worker instances, follow the directions
above for connecting workers to the master.

### How it works

Cloudbiolinux executes a sequence of installation targets listed in `fabfile.py`. Specifically:

Target `packages` installs system-wide dependencies, defined in `main.yaml` in this directory.

Target `libraries` installs packages for different languages. We only use it for python packages in `python-libs.yaml`.

Target `post_install` invokes `installer.py`, which defines installation procedures for millstone. The script clones the millstone repository from Github, sets up jbrowse, and writes a one-time bootstrap script `bootstrap.sh` located in clone directory set by `INSTALLATION_PATH`.

`bootstrap.sh` will not run until first reboot, and will run on startup after an instance is created from AMI. It generates a new password for PostgreSQL and RabbitMQ if running on master VM. This password is saved to `password.txt` along with `bootstrap.sh`.

To reconfigure, remove `BOOTSTRAPPED` file and re-run `bootstrap.sh`. Note this will generate a new password.

### Thanks

Big thanks to Changping Chen (https://github.com/ccp0101) for designing and implementing the Millstone integration with cloudbiolinux.
