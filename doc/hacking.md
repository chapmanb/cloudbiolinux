= Hacking BioLinux tips and tricks

The BioLinux tools allow building a full environment for Bioinformatics. The
design allows for flexible targets (Editions) and specializations (Flavors).

== Start with the Minimal edition

The Minimal edition is the smallest common denominator of all Editions, as it
installs the minimum of packages to bootstrap a full install. Minimal is invoked by

          fab -f $source/fabfile.py -H target_hostname -c $source/contrib/minimal/fabricrc_debian.txt install_bare:packagelist=$source/contrib/minimal/main.yaml

The main.yaml file ascertains the major editors are included, as well remote
access, version control, and the basic build system (gcc and friends). Note the
Minimal edition overwrites the (apt) sources file to make sure there are no
conflicts with user settings.

== Adding install packages

To expand on the package list you can define your own main.yaml, and pass that
in. In your main.yaml file add the meta-packages listed in
config/packages.yaml. Invoke your new package list with

          fab -f $source/fabfile.py -H target_hostname -c $source/contrib/minimal/fabricrc_debian.txt install_bare:packagelist=myproject/main.yaml

It is that simple!

If packages.yaml is not complete, you may suggest changing its contents in the
main repository. The alternative is to create your own flavor, which we will do
in a minute. The same strategy holds for the other definitions in the ./config
directory, such as for Ruby gems, Python eggs, Perl CPAN, R-CRAN etc.

== Define a Flavor

For a cross language Bio* project performance test I needed to create a special
version of BioLinux that would pull in a list of scripts and some additional
packages.  Starting from an existing edition (in this case the Minimum edition,
but it also works on top of BioNode and BioLinux editions), I created a new
flavor in ./contrib/flavor/pjotrp/biotest/biotestflavor.py, named BioTestFlavor
(note you also need an empty __init__.py file).  The flavor comes with a new
fabricrc.txt file, and a new main.yaml file.  So kicking it into submission
would look like:

          fab -f $source/fabfile.py -H target_hostname -c $source/contrib/flavor/pjotrp/biotest/fabricrc_debian.txt install_bare:packagelist=$source/contrib/flavor/pjotrp/biotest/main.yaml

The flavor module sets env.flavor (this can only happen once). For examples
see the files in ./contrib/flavor

= Tips and tricks

== Tip for checking BioLinux installation effects

To see the what a BioLinux install does to your system, store the settings of
the original (untouched) state of a VM:

1. Make a dump of the current installed package list
2. Store the /etc tree - one way is to use git in /etc

After running BioLinux you can see what has been done to your system by diffing
against the package list, and checking /etc.
