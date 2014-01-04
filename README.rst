CloudBioLinux is a build and deployment system which installs a large
selection of Bioinformatics and machine learning libraries on a bare
virtual machine (VM) image, freshly installed PC, or in the Cloud. By
default CloudBioLinux includes a large suite of tools installed through
the default distribution installer, native installers, and libraries for
Perl, R, Python, Java and Ruby.

CloudBioLinux included software packages are fully customizable. In
addition to the default configuration, we support custom configuration
builds through flavors. Flavors support overriding default package
installations, making it simple to create derived installs for specific
purposes.

CloudBioLinux is a single install route both for desktop VMs such as
`VirtualBox <http://digitizor.com/2011/01/07/virtualbox-4-0-install-ubuntu/>`_,
cloud providers such as `Amazon EC2 <http://aws.amazon.com/ec2/>`_ or
desktop machines. This works equally well for other virtual machines and
private cloud environments, including `XEN <http://xen.org/>`_, Linux
`KVM <http://www.linux-kvm.org/>`_,
`Eucalyptus <http://open.eucalyptus.com/>`_ and
`Openstack <http://www.openstack.org/>`_.

Using pre-built images
======================

Amazon
------

See the 'Getting Started with CloudBioLinux' guide on the `CloudBioLinux
website <http://cloudbiolinux.org/>`_ for a detailed description. The
short version for users familiar with Amazon is:

-  Login to the `Amazon EC2
   console <https://console.aws.amazon.com/ec2/home>`_.
-  Click Launch Instance, and choose the latest CloudBioLinux AMI from
   the `website <http://cloudbiolinux.org/>`_ in the community AMI
   section (search for 'CloudBioLinux').
-  After launching the instance, find the host details of your running
   instance from the Instances section.
-  Connect to your machine via ssh or VNC (using the Amazon PEM keys)

Installing CloudBioLinux on a local machine
===========================================

The install process for CloudBioLinux is fully automated through a
`Fabric build file <http://fabfile.org/>`_ written in Python. The Fabric
build files are useful for automating installation of scientific
software on local systems as well as Amazon cloud servers. Everything is
fully configurable through plain text YAML configuration files, and
custom build targets allow installation of a subset of the total
available packages.

Installation
------------

Retrieve the CloudBioLinux code base and install libraries and
dependencies:

::

        git clone git://github.com/chapmanb/cloudbiolinux.git
        cd cloudbiolinux
        python setup.py build
        sudo python setup.py install

Usage
-----

The basic usage specifies the hostname of a machine accessible via ssh:

::

      fab -f fabfile.py -H localhost install_biolinux

Fabric contains some other useful commandline arguments for customizing
this to your environments:

-  ``-c your_fabricrc.txt`` -- Specify the path to a fabricrc
   configuration files. This allows customization of install directories
   and other server specific details. See the default
   ``config/fabricrc.txt`` for a full list of options.

-  ``-u username`` -- The username on the remote machine, overriding the
   default of your current username.

Customization with flavors
--------------------------

CloudBioLinux normally creates a full system for bioinformatics, but can
be easily configured to install only a subset of tools through flavors:

::

      fab -f fabfile.py -H localhost install_biolinux:flavor=my_flavor

``my_flavor`` can be the name of an existing flavor in
``contrib/flavor`` or the path to a directory with customization
information. The files in your flavor directory replace those in the
standard ``config`` directory, allowing replacement of any of the
configuration files like ``main.yaml`` with customized copies.

If you desire even more control, flavors allow custom python hooks. See
``doc/hacking.md`` for more details.

Specific install targets
------------------------

You can substitute ``install_biolinux`` with more specific targets to
only build portions of CloudBioLinux:

-  ``install_biolinux:packages`` -- Install all of the defined system
   packages.
-  ``install_biolinux:libraries`` -- Install all libraries for various
   programming languages.
-  ``install_libraries:language`` -- Install libraries for a specific
   language.
-  ``install_biolinux:custom`` -- Install all custom programs.
-  ``install_custom:a_package_name`` -- Install a specific custom
   program.

Specific package installation
-----------------------------

The custom directory contains installation instructions for programs
that are not available from standard package repositories. These
instructions are written in Python using the
`Fabric <http://fabfile.org/>`_ remote deployment tool and can also be
used for installing individual packages locally on your machine. To do
this, run:

::

      fab -f fabfile.py -H localhost install_custom:your_package_name

To build and install ``your_package_name`` on the local machine. We
welcome additional custom bioinformatics package definitions for
inclusion in CloudBioLinux. ``custom/shared.py`` contains a number of
higher level functions which make it easier to write installation
instructions.

Biological data
---------------

We manage a repository of useful public biological data on an `Amazon S3
bucket <http://s3.amazonaws.com/biodata>`_. Currently this includes
whole genomes pre-indexed for a number of popular aligners. Downloading
and installing these saves a ton of time over running the indexing steps
yourself, and eases running next-generation analyses on cloud machines.

A Fabric build script is provided to install this data on your local
machine. A `biodata configuration file in YAML
format <https://github.com/chapmanb/cloudbiolinux/blob/master/config/biodata.yaml>`_,
``config/biodata.yaml``, specifies the genomes of interest and the
aligner indexes to use. The ``config/fabricrc.txt`` file specifies
details about the system and where to install the data.

The basic commandline is:

::

    fab -f data_fabfile.py -H your_machine install_data_s3

and you can pass in custom biodata and fabricrc files with:

::

    fab -f data_fabfile.py -H your_machine -c your_fabricrc.txt install_data_s3:your_biodata.yaml

In addition to downloading and preparing the data, the script will
integrate these files with a Galaxy instance by updating appropriate
Galaxy configuration files. This makes it useful for installing data to
a local or
`cloud-based <https://bitbucket.org/galaxy/galaxy-central/wiki/cloud>`_
Galaxy server.

Not all of the genomes are hosted on the S3 bucket, but are still supported. If your
genome fails to install with install_data_s3, you might be able to download the genome
from from Ensembl, etc and prepare it::


    fab -f data_fabfile.py -H your_machine -c your_fabricrc.txt install_data:your_biodata.yaml

Supported virtual environments
==============================

Vagrant VirtualBox
------------------

Vagrant allows easy deploying and connecting to VirtualBox images. The
setup is ideal for runnig CloudBioLinux on a desktop computer. Install
`VirtualBox
4.0 <http://digitizor.com/2011/01/07/virtualbox-4-0-install-ubuntu/>`_
and `vagrant <http://vagrantup.com/>`_. Then add a pre-built CloudLinux
VirtualBox images and start it up:

::

        vagrant box add biolinux_$VERSION https://s3.amazonaws.com/cloudbiolinux/biolinux_$VERSION.box
        mkdir tmp/biolinux
        cd tmp/biolinux
        vagrant init biolinux_version

(note with vagrant you need disk space - at least 3x the image size).
The created ./Vagrantfile can be edited to get a full GUI, extra RAM,
and a local IP address. Next, fire up the image with

::

        vagrant up

Once you have a running virtual machine with CloudBioLinux, connect to
it with:

::

        vagrant ssh

no passwords needed! Get root with

::

        sudo bash

Through Vagrant additional facilities are available, such as a shared
network drive. It is also possible to tweak the image (e.g. RAM/CPU
settings, and getting the all important guest additions) by firing up
virtualbox itself. For more information, see the BioLinux `Vagrant
documentation <https://github.com/chapmanb/cloudbiolinux>`_, as well as
the documentation on the `Vagrant website <http://vagrantup.com/>`_.

Amazon
------

A bare Linux image launched in Amazon EC2 is configured from another
machine, i.e. your local desktop, using ssh and cloudbiolinux. See the
Installation section for installing CloudBioLinux with fabric.

Any cloudbiolinux distribution can be used, including Ubuntu, Debian
Linux and CentOS.

1. Go to the cloudbiolinux source and edit the ``config/fabricrc.txt``,
   to match the system you plan to install on. Specifically,
   ``distribution`` and ``dist_name`` parameters specify details about
   the type of target.

2. Start an Amazon EC2 base instance and retrieve it's DNS hostname:

-  `Alestic Ubuntu images <http://alestic.com/>`_
-  `Camptocamp Debian
   images <http://www.camptocamp.com/en/infrastructure-solutions/amazon-images>`_

3. From your local machine, have CloudBioLinux install your Amazon
   instance:

   ::

       fab -f fabfile.py -H hostname -u username -i private_key_file install_biolinux

4. When finished, use the `Amazon
   console <https://console.aws.amazon.com/ec2/home>`_ to create an AMI.
   Thereafter make it public so it can be used by others.

Virtualbox
----------

See `the VirtualBox and Vagrant
documentation <https://github.com/chapmanb/cloudbiolinux/blob/master/doc/virtualbox.md>`_
for details on creating a local virtual machine from scratch with
CloudBioLinux.

OpenStack/XEN/KVM/Eucalyptus private Cloud
------------------------------------------

As long as there is an 'ssh' entry to an running VM, CloudBioLinux can
install itself.

For more on private Cloud and CloudBioLinux see ./doc/private\_cloud.md.

EC2 quickstart
==============

This provides a quick cheat sheet of commands for getting up and running
on EC2 using Amazon's command line tools.

Initial set up
--------------

The first time using EC2, you'll need to install the toolkit and
credentials for connecting on your local machine, following the `getting
started
guide <http://docs.amazonwebservices.com/AWSEC2/latest/GettingStartedGuide/>`_.

Login to your `Amazon EC2 account <http://aws.amazon.com/account/>`_ and
go to Security Credentials/X.509. Create a new certificate and download
the public ``cert-*.pem`` and ``private pk-*.pem`` files. Put these in
``~.ec2``.

Install the `ec2 api
tools <http://developer.amazonwebservices.com/connect/entry.jspa?externalID=351&categoryID=88>`_,
which require java.

Set up .zshrc/.bashrc:

::

       export EC2_PRIVATE_KEY=~/.ec2/pk-UBH43XTAWVNQMIZRAV3RP5IIBAPBIFVP.pem
       export EC2_CERT=~/.ec2/cert-UBH43XTAWVNQMIZRAV3RP5IIBAPBIFVP.pem
       export AWS_ACCESS_KEY_ID=<your access key>
       export AWS_SECRET_ACCESS_KEY=<your secret access key>

To test, you should be able to run the command:

::

       % ec2-describe-regions

Now generate a privatekey for logging in:

::

       % ec2-add-keypair yourmachine-keypair

This will produce an RSA private key. You should copy and paste this to
your .ec2 directory for future use:

::

       % vim ~/.ec2/id-yourmachine.keypair
       % chmod 600 ~/.ec2/id-yourmachine.keypair

Allow ssh and web access to your instances:

::

       % ec2-authorize default -p 22
       % ec2-authorize default -p 80

Starting an instance
--------------------

Each time you'd like to use EC2, you need to create a remote instance to
work with; the `AWS console <http://alestic.com/>`_ is useful for
managing this process.

When building from scratch with Alestic images, you will need to
increase the size of the root filesystem to fit all of the CloudBioLinux
data and libraries. This is done by starting the instance from the
commandline with:

::

       % ec2-run-instances ami-1aad5273 -k kunkel-keypair -t m1.large
                           -b /dev/sda1=:20
       % ec2-describe-instances i-0ca39764

On Ubuntu 10.04, you then need to ssh into the instance and resize the
filesystem with:

::

       % sudo resize2fs /dev/sda1

On 11.04 the resize happens automatically and this is not required.

Testing
=======

BioLinux comes with an integration testing frame work - currently based
on Vagrant. Try:

::

        cd test
        ./testing_vagrant --help

Target VMs can be listed with

::

        ./testing_vagrant --list

Build a minimal VM

::

        ./testing_vagrant Minimal

Documentation
=============

Additional documentation can be found in the `./doc
directory <https://github.com/chapmanb/cloudbiolinux>`_ in the BioLinux
source tree.

LICENSE
=======

The code is freely available under the `MIT
license <http://www.opensource.org/licenses/mit-license.html>`_.
