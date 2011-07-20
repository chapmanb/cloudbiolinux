CloudBioLinux is a build and deployment system which installs a large
selection of Bioinformatics and machine learning libraries on a bare
virtual machine (VM) image, freshly installed PC, or in the cloud. By
default CloudBioLinux includes a large suite of tools and libraries,
largely pulled from the package management system provided by the
image. In addition to the default configuration, other software
configurations are supported through, so called, 'editions' and
'flavors'. An edition is a named base install, e.g. CloudBioLinux for
Ubuntu. Such an Edition reflects shared installation properties
between different projects. Flavors support overriding base packages
for specific installations. With a Flavor a derivative of a base
Edition can be designed. For example a single web server cound be
implemented as a Flavor on top of the 'Minimal' Debian edition.
CloudBioLinux installs packages through multiple mechanisms, including
the default distribution installer, native installers, and libraries
for Perl, R, Python, JAVA and Ruby, as well as installers for special
data resources.

CloudBioLinux is designed as a single install route for both VMs on
the desktop, such as [VirtualBox][v2], and cloud providers, such as
[Amazon EC2][0], where you start with a bare bones system and
bootstrap a running instance. CloudBioLinux included software packages
are fully customizable, and different flavours of CloudBioLinux can be
configured.

CloudBioLinux provides a [Fabric build file][3], written in Python.
There is little need to understand Python, though, as most configuration
is done through configuration files. Other scripting languages can be
added too. Package selection is through YAML files in the ./config
directory.

# Using an instance

## Amazon

See the 'Getting Started with CloudBioLinux' guide on the
[CloudBioLinux website][1] for a detailed description. The short
version for users familiar with Amazon is:

* Login to the [Amazon EC2 console][2].
* Click Launch Instance, and choose the latest CloudBioLinux AMI from
  the [website][1] in the community AMI section (search for
  'biolinux').
* After launching the instance, find the host details of
  your running instance from the Instances section.
* Connect to your machine via ssh or VNC (using the Amazon PEM keys)

## Biological data

We manage a repository of useful public biological data on an
[Amazon S3 bucket][bd1]. Currently this includes whole genomes
pre-indexed for a number of popular aligners. Downloading and
installing these saves a ton of time over running the indexing steps
yourself, and eases running next-generation analyses on cloud
machines.

A Fabric build script is provided to install this data on your
local machine. A [biodata configuration file in YAML format][bd2],
`config/biodata.yaml`, specifies the genomes of interest and the
aligner indexes to use. The `config/fabricrc.txt` file specifies
details about the system and where to install the data.

The basic commandline is:

    fab -f data_fabfile.py -H your_machine install_data_s3

and you can pass in custom biodata and fabricrc files with:

    fab -f data_fabfile.py -H your_machine -c your_fabricrc.txt install_data_s3:your_biodata.yaml

In addition to downloading and preparing the data, the script will
integrate these files with a Galaxy instance by updating appropriate
Galaxy configuration files. This makes it useful for installing data
to a local or [cloud-based][bd3] Galaxy server.

[bd1]: http://s3.amazonaws.com/biodata
[bd2]: https://github.com/chapmanb/cloudbiolinux/blob/master/config/biodata.yaml
[bd3]: https://bitbucket.org/galaxy/galaxy-central/wiki/cloud

## VirtualBox with vagrant

Vagrant allows easy deploying and connecting to VirtualBox images.
Install [VirtualBox 4.0][v2] and [vagrant][v1]. Then add the pre-built
CloudLinux VirtualBox images and start it up:

        vagrant box add biolinux_version https://s3.amazonaws.com/chapmanb/biolinux_version.box
        mkdir tmp/biolinux
        cd tmp/biolinux
        vagrant init biolinux_version

(note with vagrant you need disk space - at least 3x the image size).
The created ./Vagrantfile can be edited to get a full GUI, extra RAM, and
a local IP address. Next, fire up the image with

        vagrant up

Once you have a running virtual machine with CloudBioLinux,
connect to it with:

        vagrant ssh

no passwords needed! Get root with

        sudo bash

Through Vagrant additional facilities are available, such as a shared
network drive.  It is also possible to tweak the image (e.g. RAM/CPU
settings, and getting the all important guest additions) by firing up
virtualbox itself. For more information, see the BioLinux 
[Vagrant documentation][doc], as well as the 
documentation on the [Vagrant website][v1].

# Building an image from scratch using CloudBioLinux

## Amazon

Basically a bare Linux image is configured from another machine, e.g.
your local desktop, using ssh and Fabric tools, after launching the
image in the Cloud. Multiple distributions are supported,
including Ubuntu, Debian Linux and CentOS.

1. On your local machine, install [Fabric][3]:

        sudo apt-get install python-setuptools
        sudo easy_install fabric pyyaml

2. Retrieve the cloudbiolinux code base:

        git clone git://github.com/chapmanb/cloudbiolinux.git
        cd cloudbiolinux

3. Edit the `config/fabricrc.txt` to match the system you plan to
   install on. Specifically, `distribution` and `dist_name` parameters
   specify details about the type of machine.

4. Start an Amazon EC2 base instance and get it's hostname:

   - [Alestic Ubuntu images][4]
   - [Camptocamp Debian images][4b]

5. From your local machine, start installing CloudBioLinux on your
   Amazon instance:

        fab -f fabfile.py -H hostname -u username -i private_key_file install_biolinux

6. When finished, use the [Amazon console][2] to create an AMI. Make
   it public to share it with others.

## VirtualBox with vagrant

Add a base image and boot it up; community Vagrant boxes are available from
[http://vagrantbox.es][v3]:

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

[0]: http://aws.amazon.com/ec2/
[1]: http://cloudbiolinux.org/
[2]: https://console.aws.amazon.com/ec2/home
[3]: http://fabfile.org/
[4]: http://alestic.com/
[4b]: http://www.camptocamp.com/en/infrastructure-solutions/amazon-images
[v1]: http://vagrantup.com/
[v2]: http://digitizor.com/2011/01/07/virtualbox-4-0-install-ubuntu/
[v3]: http://vagrantbox.es/

# Technical details for using build scripts

## Targets for local builds

The Fabric build files are useful for automating installation of
scientific software on local systems as well as Amazon cloud
servers. There are a number of build targets that can be used to
install a subset of the total available packages.

The main build command is used as described above, but the host
can point to a local machine or any server on your network:

      fab -f fabfile.py -H localhost -c config/fabricrc.txt install_biolinux

The `config/fabricrc.txt` configuration file specifies install
directories and other server specific details.

With this basic command, you can substitute `install_biolinux` with
serveral more specific targets:

* `install_biolinux:packages` -- Install all of the defined system
  packages.
* `install_biolinux:libraries` -- Install all libraries for various
  programming languages.
* `install_libraries:language` -- Install libraries for a specific
  language.
* `install_biolinux:custom` -- Install all custom programs.
* `install_custom:your_package_name` -- Install a specific custom
   program.

### Custom package installs

The custom directory contains installation instructions for programs that are
not available from standard package repositories. These instructions are written
in Python using the [Fabric][3] remote deployment tool and can also be used for
installing individual packages locally on your machine. To do this, run:

      fab -f fabfile.py -H localhost install_custom:your_package_name

To build and install `your_package_name` on the local machine. We welcome
additional custom bioinformatics package definitions for inclusion in
CloudBioLinux. `custom/shared.py` contains a number of higher level functions
which make it easier to write installation instructions.

### Using flavors

A Flavor is a specialization of an installation edition. A flavor can filter on
packages and add packages and other software. To use a flavor, pass it
in on the command line, e.g.

      fab -f fabfile.py -H localhost install_biolinux:flavor=my_flavor

And, like with the other options, a flavor can also be defined in the
config file. A command line parameter, however, is preferred.

### Rolling your own

BioLinux normally creates a full system for Bioinformatics. The tools
provided here to create such an environment are rather generic. So, if you want to
install packages on any system, be it desktop, server or VM, the
BioLinux tool set can be used to role your own. The first step it to us
install_biolinux with a named package list, e.g.

      fab -f fabfile.py -H localhost install_biolinux:packagelist=./contrib/minimal/main.yaml

where install_biolinux is the fab entry point to a non-specific install,
starting from the package list in the passed in packagelist. For more
information see the section 'Rolling your own' in the BioLinux vagrant documentation in [./doc/vagrant.md][doc].

## EC2 quickstart

This provides a quick cheat sheet of commands for getting up and running on EC2 using
Amazon's command line tools.

### Initial set up

The first time using EC2, you'll need to install the toolkit and credentials
for connecting on your local machine, following the [getting started guide][qs1].

Login to your [Amazon EC2 account][qs2] and go to Security Credentials/X.509.
Create a new certificate and download the public `cert-*.pem` and
`private pk-*.pem` files. Put these in `~.ec2`.

Install the [ec2 api tools][qs3], which require java.

Set up .zshrc/.bashrc:

       export EC2_PRIVATE_KEY=~/.ec2/pk-UBH43XTAWVNQMIZRAV3RP5IIBAPBIFVP.pem
       export EC2_CERT=~/.ec2/cert-UBH43XTAWVNQMIZRAV3RP5IIBAPBIFVP.pem
       export AWS_ACCESS_KEY_ID=<your access key>
       export AWS_SECRET_ACCESS_KEY=<your secret access key>

To test, you should be able to run the command:

       % ec2-describe-regions

Now generate a privatekey for logging in:

       % ec2-add-keypair yourmachine-keypair

This will produce an RSA private key. You should copy and paste this to your
.ec2 directory for future use:

       % vim ~/.ec2/id-yourmachine.keypair
       % chmod 600 ~/.ec2/id-yourmachine.keypair

Allow ssh and web access to your instances:

       % ec2-authorize default -p 22
       % ec2-authorize default -p 80

[qs1]: http://docs.amazonwebservices.com/AWSEC2/latest/GettingStartedGuide/
[qs2]: http://aws.amazon.com/account/
[qs3]: http://developer.amazonwebservices.com/connect/entry.jspa?externalID=351&categoryID=88

### Starting an instance

Each time you'd like to use EC2, you need to create a remote instance to work
with; the [AWS console][4] is useful for managing this process.

When building from scratch with Alestic images, you will need to
increase the size of the root filesystem to fit all of the
CloudBioLinux data and libraries. This is done by starting the
instance from the commandline with:

       % ec2-run-instances ami-1aad5273 -k kunkel-keypair -t m1.large
                           -b /dev/sda1=:20
       % ec2-describe-instances i-0ca39764

On Ubuntu 10.04, you then need to ssh into the instance and resize the
filesystem with:

       % sudo resize2fs /dev/sda1

On 11.04 the resize happens automatically and this is not required.

# Testing

BioLinux comes with an integration testing frame work - currently
based on Vagrant. Try:

        cd test
        ./testing_vagrant --help

# Documentation

Additional documentation can be found in the [./doc directory][doc] in the
BioLinux source tree.

[doc]: https://github.com/chapmanb/cloudbiolinux

# LICENSE

The code is freely available under the [MIT license][l1].

[l1]: http://www.opensource.org/licenses/mit-license.html
