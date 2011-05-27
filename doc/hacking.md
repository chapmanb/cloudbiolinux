= Hacking BioLinux tips and tricks

The BioLinux tools allow building a full environment for Bioinformatics. The
design allows for flexible targets (Editions) and specializations (Flavors).

== The Minimal edition

The Minimal edition is the smallest common denominator of all Editions, as it
installs the minimum of packages to bootstrap a full install. Minimal is invoked by

          fab -f $source/fabfile.py -H target_hostname -c $source/contrib/minimal/fabricrc_debian.txt install_bare:packagelist=$source/contrib/minimal/main.yaml

The main.yaml file ascertains the major editors are included, as well remote
access, version control, and the basic build system (gcc and friends). Note the
Minimal edition overwrites the (apt) sources file to make sure there are no
conflicts with user settings.

= Tips and tricks

== Tip for checking BioLinux installation effects

To see the what a BioLinux install does to your system, store the settings of
the original (untouched) state of a VM:

1. Make a dump of the current installed package list
2. Store the /etc tree - one way is to use git in /etc

After running BioLinux you can see what has been done to your system by diffing
against the package list, and checking /etc.
